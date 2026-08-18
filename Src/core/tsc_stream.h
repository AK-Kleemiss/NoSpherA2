#pragma once
//
// Streaming, overlapped writer for scattering-factor tables.
//
// Two problems, one mechanism:
//
//   memory - the structure-factor array is scatterers x reflections x 16 bytes,
//            which is 21 GB for a 8500-atom protein at 157k reflections and is
//            by far the largest allocation in a run. The tscb payload is
//            reflection-major, so blocks of reflections can be emitted in order
//            and only a few blocks ever need to be resident.
//
//   time   - writing gigabytes is not free. A dedicated writer thread lets the
//            next block be computed while the previous one goes to disk.
//
// The producer side is deliberately a plain callback, so every kind of
// scattering factor - SALTED, spherical Thakkar, wavefunction-based, and the
// electron-diffraction conversion - can use the same path rather than each
// growing its own buffering.
//
#include "pch.h"
#include "tsc_block.h"

#include <condition_variable>
#include <mutex>
#include <thread>
#include <map>

template <typename numtype_index, typename numtype>
class tsc_stream_writer
{
public:
    using block_type = tsc_block<numtype_index, numtype>;
    using index_block = std::vector<std::vector<numtype_index>>;

    // queue_depth bounds the memory: at most this many blocks are resident
    // beyond the one being filled.
    tsc_stream_writer(const std::filesystem::path& name,
                      const block_type& prototype,
                      const std::size_t reflection_count,
                      const std::size_t queue_depth = 2)
        : out_(name, std::ios::binary | std::ios::trunc),
          queue_depth_(queue_depth == 0 ? 1 : queue_depth)
    {
        if (!out_)
            throw std::runtime_error("Failed to open TSCB file for streaming");
        prototype.write_tscb_prologue(out_, reflection_count);
        worker_ = std::thread([this] { drain(); });
    }

    // The production entry point: only the scatterer list is needed, so no
    // full-sized block is ever built.
    tsc_stream_writer(const std::filesystem::path& name,
                      const ScattererLabels& scatterers,
                      const std::string& header,
                      const std::size_t reflection_count,
                      const std::size_t queue_depth = 2)
        : out_(name, std::ios::binary | std::ios::trunc),
          queue_depth_(queue_depth == 0 ? 1 : queue_depth)
    {
        if (!out_)
            throw std::runtime_error("Failed to open TSCB file for streaming");
        block_type::write_tscb_prologue(out_, scatterers, header, reflection_count);
        worker_ = std::thread([this] { drain(); });
    }

    ~tsc_stream_writer() { try { finish(); } catch (...) { } }

    tsc_stream_writer(const tsc_stream_writer&) = delete;
    tsc_stream_writer& operator=(const tsc_stream_writer&) = delete;

    // Blocks may be produced out of order; they are written in order of id.
    void submit(const std::size_t id, index_block idx, cvec2 sf_block)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        space_.wait(lock, [this] { return pending_.size() < queue_depth_ || failed_; });
        rethrow_if_failed();
        pending_.emplace(id, item{std::move(idx), std::move(sf_block)});
        lock.unlock();
        work_.notify_one();
    }

    void finish()
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (done_) return;
            done_ = true;
        }
        work_.notify_all();
        if (worker_.joinable()) worker_.join();
        std::lock_guard<std::mutex> lock(mutex_);
        rethrow_if_failed();
        out_.flush();
    }

private:
    struct item { index_block idx; cvec2 sf; };

    void rethrow_if_failed() const
    {
        if (error_) std::rethrow_exception(error_);
    }

    void drain()
    {
        try
        {
            for (;;)
            {
                item next;
                {
                    std::unique_lock<std::mutex> lock(mutex_);
                    work_.wait(lock, [this] {
                        return done_ || pending_.count(next_id_) != 0;
                    });
                    const auto found = pending_.find(next_id_);
                    if (found == pending_.end())
                    {
                        // Nothing more can arrive for this id once finish() ran,
                        // and anything left is out of order, which is a caller bug.
                        if (done_ && pending_.empty()) return;
                        if (done_) throw std::runtime_error(
                            "tsc_stream_writer: missing reflection block " + std::to_string(next_id_));
                        continue;
                    }
                    next = std::move(found->second);
                    pending_.erase(found);
                    ++next_id_;
                }
                space_.notify_one();
                block_type::write_tscb_reflection_block(out_, next.idx, next.sf);
            }
        }
        catch (...)
        {
            std::lock_guard<std::mutex> lock(mutex_);
            error_ = std::current_exception();
            failed_ = true;
            space_.notify_all();
        }
    }

    std::ofstream out_;
    std::size_t queue_depth_;
    std::map<std::size_t, item> pending_;
    std::size_t next_id_ = 0;
    bool done_ = false;
    bool failed_ = false;
    std::exception_ptr error_;
    std::mutex mutex_;
    std::condition_variable work_, space_;
    std::thread worker_;
};
