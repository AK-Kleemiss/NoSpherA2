#pragma once
#include "convenience.h"
#include <cassert>

// The XCW I tensor is nr blocks of complex values, one per reflection, holding only the
// (mu, nu) pairs that survived the overlap screening - kept of the nmo*(nmo+1)/2 - in the
// order the pair list in the header gives. The header also says whether the elements are
// single or double: a tensor built in single precision is stored as it was computed. Every
// consumer walks it in reflection order, so only a window need be resident:
//     resident = window * kept * (8 or 16) bytes.
// C stdio, not ifstream, for read speed; only the 64-bit seek needs an #if.
class i_tensor_file
{
public:
    static size_t block_bytes(const size_t kept, const bool single)
    {
        return kept * (single ? sizeof(std::complex<float>) : sizeof(cdouble));
    }

    // size_t, not int: the element count wraps at 2^31 elements (34 GB) and truncates.
    static size_t total_bytes(const int nr, const size_t kept, const bool single)
    {
        return static_cast<size_t>(nr) * block_bytes(kept, single);
    }

    // Whether p holds a complete tensor for this problem; kept and single receive its pair
    // count and element type. A file in the previous layout, every pair stored in double, is
    // refused with a note: its blocks cannot be read into the compact one.
    static bool matches(const std::filesystem::path &p, const int nr, const int nmo, size_t &kept, bool &single)
    {
        kept = 0;
        single = false;
        FILE *f = fopen(p.string().c_str(), "rb");
        if (!f) return false;
        int64_t h[5] = { 0, 0, 0, 0, 0 };
        const bool got = fread(h, sizeof(int64_t), 5, f) == 5;
        fclose(f);
        if (!got) return false;
        if (h[0] == magic_full_)
        {
            std::cout << "I tensor " << p.string() << " stores every pair, the previous layout; it is rebuilt" << std::endl;
            return false;
        }
        if (h[0] != magic_ || h[1] != nr || h[2] != nmo || h[3] <= 0) return false;
        if (h[4] != static_cast<int64_t>(sizeof(cdouble)) && h[4] != static_cast<int64_t>(sizeof(std::complex<float>))) return false;
        const bool s = h[4] == static_cast<int64_t>(sizeof(std::complex<float>));
        std::error_code ec;
        const auto on_disk = std::filesystem::file_size(p, ec);
        if (ec || static_cast<size_t>(on_disk) < header_bytes(static_cast<size_t>(h[3])) + total_bytes(nr, static_cast<size_t>(h[3]), s))
            return false;
        kept = static_cast<size_t>(h[3]);
        single = s;
        return true;
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
    void create(const std::filesystem::path &p, const int nr, const int nmo, const ivec &pair_mu, const ivec &pair_nu, const bool single)
    {
        close();
        nr_ = nr; nmo_ = nmo; kept_ = pair_mu.size(); single_ = single;
        pair_mu_ = pair_mu; pair_nu_ = pair_nu;
        f_ = fopen(p.string().c_str(), "wb+");
        if (!f_) throw std::runtime_error("i_tensor_file: cannot create " + p.string());
        const int64_t h[5] = { magic_, nr, nmo, static_cast<int64_t>(kept_), static_cast<int64_t>(elem_bytes()) };
        fwrite(h, sizeof(int64_t), 5, f_);
        fwrite(pair_mu_.data(), sizeof(int), kept_, f_);
        fwrite(pair_nu_.data(), sizeof(int), kept_, f_);
        path_ = p;
    }

    // Either element type in, the file's type out. Callers serialise these; the conversion
    // scratch is one member.
    void write_block(const int r, const cdouble *block)
    {
        if (!single_) { write_raw(r, block); return; }
        scratch32_.resize(kept_);
        for (size_t i = 0; i < kept_; i++)
            scratch32_[i] = std::complex<float>(static_cast<float>(block[i].real()), static_cast<float>(block[i].imag()));
        write_raw(r, scratch32_.data());
    }
    void write_block(const int r, const std::complex<float> *block)
    {
        if (single_) { write_raw(r, block); return; }
        scratch64_.resize(kept_);
        for (size_t i = 0; i < kept_; i++) scratch64_[i] = cdouble(block[i].real(), block[i].imag());
        write_raw(r, scratch64_.data());
    }

    void finish_write() { if (f_) fflush(f_); }

    void open(const std::filesystem::path &p, const size_t window_blocks)
    {
        close();
        f_ = fopen(p.string().c_str(), "rb");
        if (!f_) throw std::runtime_error("i_tensor_file: cannot open " + p.string());
        int64_t h[5] = { 0, 0, 0, 0, 0 };
        if (fread(h, sizeof(int64_t), 5, f_) != 5 || h[0] != magic_)
            throw std::runtime_error("i_tensor_file: " + p.string() + " is not an I tensor in the compact layout");
        if (h[4] != static_cast<int64_t>(sizeof(cdouble)) && h[4] != static_cast<int64_t>(sizeof(std::complex<float>)))
            throw std::runtime_error("i_tensor_file: " + p.string() + " has an unknown element size");
        single_ = h[4] == static_cast<int64_t>(sizeof(std::complex<float>));
        // These header fields size an allocation and compute seek offsets, so a stale or
        // truncated file must be rejected here, not read past the end.
        if (h[1] <= 0 || h[2] <= 0 || h[3] <= 0)
            throw std::runtime_error("i_tensor_file: " + p.string() +
                " has a non-positive dimension in its header");
        if (h[1] > (1LL << 31) || h[2] > (1LL << 20) || h[3] > h[2] * (h[2] + 1) / 2)
            throw std::runtime_error("i_tensor_file: " + p.string() +
                " claims implausible dimensions (nr " + std::to_string(h[1]) +
                ", nmo " + std::to_string(h[2]) + ", pairs " + std::to_string(h[3]) + ")");
        nr_ = static_cast<int>(h[1]); nmo_ = static_cast<int>(h[2]);
        kept_ = static_cast<size_t>(h[3]);
        pair_mu_.resize(kept_); pair_nu_.resize(kept_);
        if (fread(pair_mu_.data(), sizeof(int), kept_, f_) != kept_ ||
            fread(pair_nu_.data(), sizeof(int), kept_, f_) != kept_)
            throw std::runtime_error("i_tensor_file: truncated pair list in " + p.string());
        // The pair list indexes the density matrix in the SCF walks, so a damaged one would
        // read out of bounds there rather than fail here.
        for (size_t k = 0; k < kept_; k++)
        {
            const int mu = pair_mu_[k], nu = pair_nu_[k];
            const bool ordered = k == 0 || pair_mu_[k - 1] < mu || (pair_mu_[k - 1] == mu && pair_nu_[k - 1] < nu);
            if (mu < 0 || mu > nu || nu >= nmo_ || !ordered)
                throw std::runtime_error("i_tensor_file: " + p.string() + " has a corrupt pair list at entry " + std::to_string(k));
        }
        std::error_code ec;
        const auto on_disk = std::filesystem::file_size(p, ec);
        if (!ec)
        {
            const size_t expect = header_bytes(kept_) + total_bytes(nr_, kept_, single_);
            if (static_cast<size_t>(on_disk) != expect)
                throw std::runtime_error("i_tensor_file: " + p.string() + " is " +
                    std::to_string(on_disk) + " bytes, but its header implies " +
                    std::to_string(expect) + " - truncated or written by another format");
        }
        path_ = p;
        set_window(window_blocks);
    }

    void set_window(size_t window_blocks)
    {
        if (window_blocks < 1) window_blocks = 1;
        if (window_blocks > static_cast<size_t>(nr_)) window_blocks = static_cast<size_t>(nr_);
        window_blocks_ = window_blocks;
        if (single_) window32_.assign(window_blocks_ * kept_, std::complex<float>{});
        else window_.assign(window_blocks_ * kept_, cdouble{});
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
        const size_t want = n * kept_;
        const size_t got = single_ ? fread(window32_.data(), sizeof(std::complex<float>), want, f_)
                                   : fread(window_.data(), sizeof(cdouble), want, f_);
        if (got != want)
            throw std::runtime_error("i_tensor_file: short read at reflection " + std::to_string(r0));
        loaded_first_ = r0; loaded_last_ = r1;
    }

    // assert, not throw: called inside an omp for, and an exception leaving an OpenMP
    // structured block is undefined behaviour. Real bounds checking happens in load().
    const cdouble *block(const int r) const
    {
        assert(!single_ && r >= loaded_first_ && r < loaded_last_ && "reflection outside the loaded window");
        return window_.data() + static_cast<size_t>(r - loaded_first_) * kept_;
    }
    const std::complex<float> *block32(const int r) const
    {
        assert(single_ && r >= loaded_first_ && r < loaded_last_ && "reflection outside the loaded window");
        return window32_.data() + static_cast<size_t>(r - loaded_first_) * kept_;
    }

    int nr() const { return nr_; }
    int nmo() const { return nmo_; }
    size_t kept() const { return kept_; }
    bool single() const { return single_; }
    const ivec &pair_mu() const { return pair_mu_; }
    const ivec &pair_nu() const { return pair_nu_; }
    size_t window_blocks() const { return window_blocks_; }
    const std::filesystem::path &path() const { return path_; }

private:
    static constexpr int64_t magic_ = 0x4E4132495F544E32LL;        // "NA2I_TN2"
    static constexpr int64_t magic_full_ = 0x4E4132495F54454ELL;   // "NA2I_TEN", every pair stored

    static size_t header_bytes(const size_t kept)
    {
        return 5 * sizeof(int64_t) + 2 * kept * sizeof(int);
    }

    size_t elem_bytes() const { return single_ ? sizeof(std::complex<float>) : sizeof(cdouble); }

    size_t offset_of(const int r) const
    {
        return header_bytes(kept_) + static_cast<size_t>(r) * kept_ * elem_bytes();
    }

    template <typename T>
    void write_raw(const int r, const T *block)
    {
        seek(offset_of(r));
        if (fwrite(block, sizeof(T), kept_, f_) != kept_)
            throw std::runtime_error("i_tensor_file: short write on block " + std::to_string(r));
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
    size_t kept_ = 0;
    bool single_ = false;
    ivec pair_mu_, pair_nu_;
    size_t window_blocks_ = 0;
    cvec window_, scratch64_;
    std::vector<std::complex<float>> window32_, scratch32_;
    int loaded_first_ = -1, loaded_last_ = -1;
};
