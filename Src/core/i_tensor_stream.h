#pragma once
#include "convenience.h"

// The XCW I tensor is nr_small blocks of nmo*(nmo+1)/2 complex doubles, one block
// per reflection, and every consumer walks it in reflection order: calc_F_calc
// reads block r to produce F_calc[0][r], calc_perturb reads block r again to
// accumulate an nmo x nmo matrix. Nothing ever needs two blocks at once.
//
// That is the same shape as a .tscb - reflection-major records consumed in order,
// reducing into something small - and it admits the same treatment. Keeping the
// tensor on disk and holding a window of blocks makes memory a choice rather than
// a function of the reflection count:
//
//     resident = window * nmo*(nmo+1)/2 * 16 bytes
//
// instead of nr_small times the same block. For sto-3g on the P1 test that is
// 276 MB either way and irrelevant; the point is the basis-set scaling, which is
// quadratic in nmo, so a triple-zeta basis on the same data is 25x that and a
// larger structure is out of reach entirely.
//
// C stdio, not std::ifstream: measured on a 3.5 GB file, fread reads at
// 2.85-2.90 GB/s against 1.26-1.41 for ifstream, and chunking the stream does not
// close the gap - the stream itself is the cost. Only the 64-bit seek needs an
// #if, because fseek takes a long.
class i_tensor_file
{
public:
    // Bytes of one reflection's block.
    static size_t block_bytes(const int nmo)
    {
        return static_cast<size_t>(nmo) * (nmo + 1) / 2 * sizeof(cdouble);
    }

    // What the whole tensor would cost resident. size_t, deliberately: the
    // original header wrote the element count as an int, which wraps at 2^31
    // elements - 34 GB - and silently truncates the file for exactly the
    // structures this class exists to serve.
    static size_t total_bytes(const int nr, const int nmo)
    {
        return static_cast<size_t>(nr) * block_bytes(nmo);
    }

    ~i_tensor_file() { close(); }

    void close()
    {
        if (f_) { fclose(f_); f_ = NULL; }
        window_.clear();
        window_.shrink_to_fit();
    }

    // ---- writing -------------------------------------------------------
    // Blocks arrive in whatever order the workers finish them, so the writer
    // holds them by index and seeks. The file is laid out reflection-major, so
    // block r starts at header_bytes + r * block_bytes.
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

    // ---- reading -------------------------------------------------------
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

    // Read reflections [r0, r1) into the window. Call from one thread; the
    // pointers handed out by block() stay valid until the next load.
    void load(const int r0, const int r1)
    {
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

    const cdouble *block(const int r) const
    {
        return window_.data() + static_cast<size_t>(r - loaded_first_) * packed_;
    }

    int nr() const { return nr_; }
    int nmo() const { return nmo_; }
    size_t packed() const { return packed_; }
    size_t window_blocks() const { return window_blocks_; }
    const std::filesystem::path &path() const { return path_; }

private:
    static const size_t header_bytes_ = 4 * sizeof(int64_t);

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
