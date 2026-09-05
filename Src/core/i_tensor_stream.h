#pragma once
#include "convenience.h"
#include <cassert>

// The XCW I tensor is nr blocks of nmo*(nmo+1)/2 complex doubles, one per reflection, and
// every consumer walks it in reflection order, so only a window need be resident:
//     resident = window * nmo*(nmo+1)/2 * 16 bytes,  quadratic in nmo.
// C stdio, not ifstream, for read speed; only the 64-bit seek needs an #if.
class i_tensor_file
{
public:
    static size_t block_bytes(const int nmo)
    {
        return static_cast<size_t>(nmo) * (nmo + 1) / 2 * sizeof(cdouble);
    }

    // size_t, not int: the element count wraps at 2^31 elements (34 GB) and truncates.
    static size_t total_bytes(const int nr, const int nmo)
    {
        return static_cast<size_t>(nr) * block_bytes(nmo);
    }

    i_tensor_file() = default;
    // Owns a FILE*: a copy would fclose the same handle twice and fight over the position.
    i_tensor_file(const i_tensor_file &) = delete;
    i_tensor_file &operator=(const i_tensor_file &) = delete;
    ~i_tensor_file() { close(); }

    void close()
    {
        if (f_) { fclose(f_); f_ = NULL; }
        window_.clear();
        window_.shrink_to_fit();
    }

    // Blocks arrive in whatever order the workers finish them, so the writer seeks by index.
    void create(const std::filesystem::path &p, const int nr, const int nmo)
    {
        close();
        nr_ = nr; nmo_ = nmo; packed_ = static_cast<size_t>(nmo) * (nmo + 1) / 2;
        f_ = fopen(p.string().c_str(), "wb+");
        if (!f_) throw std::runtime_error("i_tensor_file: cannot create " + p.string());
        const int64_t magic = 0x4E4132495F54454ELL;   // "NA2I_TEN"
        const int64_t nr64 = nr, nmo64 = nmo, packed64 = static_cast<int64_t>(packed_);
        fwrite(&magic, sizeof(magic), 1, f_);
        fwrite(&nr64, sizeof(nr64), 1, f_);
        fwrite(&nmo64, sizeof(nmo64), 1, f_);
        fwrite(&packed64, sizeof(packed64), 1, f_);
        path_ = p;
    }

    void write_block(const int r, const cdouble *block)
    {
        seek(offset_of(r));
        if (fwrite(block, sizeof(cdouble), packed_, f_) != packed_)
            throw std::runtime_error("i_tensor_file: short write on block " + std::to_string(r));
    }

    void finish_write() { if (f_) fflush(f_); }

    void open(const std::filesystem::path &p, const size_t window_blocks)
    {
        close();
        f_ = fopen(p.string().c_str(), "rb");
        if (!f_) throw std::runtime_error("i_tensor_file: cannot open " + p.string());
        int64_t magic = 0, nr64 = 0, nmo64 = 0, packed64 = 0;
        if (fread(&magic, sizeof(magic), 1, f_) != 1 || magic != 0x4E4132495F54454ELL)
            throw std::runtime_error("i_tensor_file: " + p.string() + " is not an I tensor");
        if (fread(&nr64, sizeof(nr64), 1, f_) != 1 ||
            fread(&nmo64, sizeof(nmo64), 1, f_) != 1 ||
            fread(&packed64, sizeof(packed64), 1, f_) != 1)
            throw std::runtime_error("i_tensor_file: truncated header");
        // These header fields size an allocation and compute seek offsets, so a stale or
        // truncated file must be rejected here, not read past the end.
        if (nr64 <= 0 || nmo64 <= 0 || packed64 <= 0)
            throw std::runtime_error("i_tensor_file: " + p.string() +
                " has a non-positive dimension in its header");
        if (nr64 > (1LL << 31) || nmo64 > (1LL << 20))
            throw std::runtime_error("i_tensor_file: " + p.string() +
                " claims implausible dimensions (nr " + std::to_string(nr64) +
                ", nmo " + std::to_string(nmo64) + ")");
        const int64_t packed_expected = nmo64 * (nmo64 + 1) / 2;
        if (packed64 != packed_expected)
            throw std::runtime_error("i_tensor_file: " + p.string() + " header is inconsistent: "
                "packed " + std::to_string(packed64) + " but nmo " + std::to_string(nmo64) +
                " implies " + std::to_string(packed_expected));
        std::error_code ec;
        const auto on_disk = std::filesystem::file_size(p, ec);
        if (!ec)
        {
            const size_t expect = header_bytes_ +
                static_cast<size_t>(nr64) * static_cast<size_t>(packed64) * sizeof(cdouble);
            if (static_cast<size_t>(on_disk) != expect)
                throw std::runtime_error("i_tensor_file: " + p.string() + " is " +
                    std::to_string(on_disk) + " bytes, but its header implies " +
                    std::to_string(expect) + " - truncated or written by another format");
        }
        nr_ = static_cast<int>(nr64); nmo_ = static_cast<int>(nmo64);
        packed_ = static_cast<size_t>(packed64);
        path_ = p;
        set_window(window_blocks);
    }

    void set_window(size_t window_blocks)
    {
        if (window_blocks < 1) window_blocks = 1;
        if (window_blocks > static_cast<size_t>(nr_)) window_blocks = static_cast<size_t>(nr_);
        window_blocks_ = window_blocks;
        window_.assign(window_blocks_ * packed_, cdouble{});
        loaded_first_ = -1; loaded_last_ = -1;
    }

    // Call from one thread; pointers handed out by block() stay valid until the next load.
    void load(const int r0, const int r1)
    {
        if (r0 < 0 || r1 < r0 || r1 > nr_)
            throw std::runtime_error("i_tensor_file: load(" + std::to_string(r0) + ", " +
                std::to_string(r1) + ") is outside 0.." + std::to_string(nr_));
        const size_t n = static_cast<size_t>(r1 - r0);
        if (n > window_blocks_)
            throw std::runtime_error("i_tensor_file: load of " + std::to_string(n) +
                                     " blocks exceeds the window");
        seek(offset_of(r0));
        const size_t want = n * packed_;
        if (fread(window_.data(), sizeof(cdouble), want, f_) != want)
            throw std::runtime_error("i_tensor_file: short read at reflection " + std::to_string(r0));
        loaded_first_ = r0; loaded_last_ = r1;
    }

    // assert, not throw: called inside an omp for, and an exception leaving an OpenMP
    // structured block is undefined behaviour. Real bounds checking happens in load().
    const cdouble *block(const int r) const
    {
        assert(r >= loaded_first_ && r < loaded_last_ && "reflection outside the loaded window");
        return window_.data() + static_cast<size_t>(r - loaded_first_) * packed_;
    }

    int nr() const { return nr_; }
    int nmo() const { return nmo_; }
    size_t packed() const { return packed_; }
    size_t window_blocks() const { return window_blocks_; }
    const std::filesystem::path &path() const { return path_; }

private:
    static constexpr size_t header_bytes_ = 4 * sizeof(int64_t);

    size_t offset_of(const int r) const
    {
        return header_bytes_ + static_cast<size_t>(r) * packed_ * sizeof(cdouble);
    }

    void seek(const size_t off)
    {
#ifdef _WIN32
        if (_fseeki64(f_, static_cast<__int64>(off), SEEK_SET) != 0)
#else
        if (fseeko(f_, static_cast<off_t>(off), SEEK_SET) != 0)
#endif
            throw std::runtime_error("i_tensor_file: seek failed");
    }

    FILE *f_ = NULL;
    std::filesystem::path path_;
    int nr_ = 0, nmo_ = 0;
    size_t packed_ = 0;
    size_t window_blocks_ = 0;
    cvec window_;
    int loaded_first_ = -1, loaded_last_ = -1;
};
