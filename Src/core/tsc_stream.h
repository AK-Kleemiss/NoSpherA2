#pragma once
// Streamed tscb writer. The file is reflection-major, so a block of reflections can be
// written and forgotten; nothing later depends on anything earlier and we never seek back.
// One writer thread in drain() is the only thread that touches out_; mutex_ guards
// pending_, next_id_, done_, failed_ and error_.
// Peak memory is queue_depth * scatterers * reflections-per-block * 16 bytes, so submit()
// blocking on a full queue is that bound, not an inefficiency.
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

    // id must run 0, 1, 2, ... with no gaps or the writer stalls; blocks may arrive in any
    // order. idx and sf_block are moved from and owned by the writer thread until finish().
    void submit(const std::size_t id, index_block idx, cvec2 sf_block)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        space_.wait(lock, [this] { return pending_.size() < queue_depth_ || failed_; });
        rethrow_if_failed();
        pending_.emplace(id, item{ std::move(idx), std::move(sf_block) });
        lock.unlock();       // release first, so the woken writer does not
        work_.notify_one();  // immediately block on a mutex we still hold
    }

    // Waits for every submitted block, joins the writer, rethrows its error; idempotent.
    void finish()
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (done_) return;
            done_ = true;    // tells the writer that no more blocks are coming
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
                    // Waiting on next_id_, not on non-empty, is what restores file order.
                    work_.wait(lock, [this] {
                        return done_ || pending_.count(next_id_) != 0;
                    });
                    const auto found = pending_.find(next_id_);
                    if (found == pending_.end())
                    {
                        if (done_ && pending_.empty())
                            return;                        // finished cleanly
                        if (done_)
                            throw std::runtime_error(      // a gap in the ids;
                                "tsc_stream_writer: missing reflection block "
                                + std::to_string(next_id_));  // waiting would hang
                        continue;                          // spurious wake-up
                    }
                    next = std::move(found->second);
                    pending_.erase(found);
                    ++next_id_;
                }                     // mutex released here, before the slow part
                space_.notify_one();  // a slot freed: a waiting producer may proceed

                block_type::write_tscb_reflection_block(out_, next.idx, next.sf);
            }
        }
        catch (...)
        {
            // An exception escaping a thread terminates the program; park it for the producer.
            std::lock_guard<std::mutex> lock(mutex_);
            error_ = std::current_exception();
            failed_ = true;
            space_.notify_all();
        }
    }

    std::ofstream out_;                     // only the writer thread touches this
    std::size_t queue_depth_;               // how many blocks may be resident
    std::map<std::size_t, item> pending_;   // ordered by id, so order is restored
    std::size_t next_id_ = 0;               // the id the writer will accept next
    bool done_ = false;                     // finish() has been called
    bool failed_ = false;                   // the writer thread hit an error
    std::exception_ptr error_;              // ... and this is what it was
    std::mutex mutex_;                      // guards everything above
    std::condition_variable work_, space_;  // "block arrived" / "slot freed"
    std::thread worker_;
};
