#pragma once
//
// ============================================================================
//  Writing a scattering-factor table in pieces, while the next piece is computed
// ============================================================================
//
// WHY THIS EXISTS
//
// A tsc table holds one complex number per (scatterer, reflection). For a
// protein that is enormous: 8,566 atoms x 293,000 reflections x 16 bytes is
// 40 GB, and the original code built all of it in memory before writing
// anything. That alone decided whether a structure could be processed.
//
// Two facts make a better approach possible:
//
//   1. The tscb file is laid out REFLECTION BY REFLECTION:
//
//          reflection 0:  h k l   sf(atom0)  sf(atom1) ... sf(atomN)
//          reflection 1:  h k l   sf(atom0)  sf(atom1) ... sf(atomN)
//          ...
//
//      so a block of reflections can be finished, written, and forgotten before
//      the next block starts. Nothing later in the file depends on anything
//      earlier, and we never seek backwards.
//
//   2. Writing gigabytes is slow, and while the disk works the CPU idles. Put
//      the writing on its own thread and the next block is computed meanwhile.
//
// ============================================================================
//  HOW IT WORKS - written for readers new to threaded code
// ============================================================================
//
// This is the classic "producer / consumer" arrangement.
//
//   The PRODUCER is ordinary calculation code. It computes a block of
//   reflections and calls submit(). It never touches the file.
//
//   The CONSUMER is one extra thread, started in the constructor, which sits in
//   drain() taking blocks and writing them. It is the ONLY thread that ever
//   touches the file - which is what keeps the file consistent without locking
//   around the writing itself.
//
//   Between them sits a QUEUE (pending_): the producer puts blocks in, the
//   consumer takes them out.
//
// Three problems come with that arrangement. Each has a specific mechanism here.
// If you read nothing else, read these.
//
// PROBLEM 1 - two threads touching the same data at once.
//
//   If the producer inserts into the queue while the consumer removes from it,
//   the queue's internals are corrupted. Not "occasionally wrong numbers" -
//   arbitrary memory damage.
//
//   Solution: a MUTEX (mutex_), a token only one thread may hold. Every piece
//   of code touching pending_, next_id_, done_ or error_ takes that token first
//   (std::lock_guard / std::unique_lock, which release it automatically at the
//   end of the enclosing block). While one thread holds it, others wait.
//
// PROBLEM 2 - waiting without burning the CPU.
//
//   The consumer has nothing to do until a block arrives. Spinning in a loop
//   asking "anything yet?" would consume a whole core doing nothing.
//
//   Solution: a CONDITION VARIABLE (work_, space_). wait() puts the thread to
//   sleep using no CPU, and releases the mutex while it sleeps. Another thread
//   calls notify_one() to wake it. The lambda passed to wait() is the condition
//   being waited for, and wait() re-checks it on every wake-up - threads can
//   wake spuriously, and another thread may have taken the work first. Never
//   replace it with a bare wait().
//
//   There are two because there are two different things to wait for:
//       work_  - the consumer waits for a block to arrive
//       space_ - the producer waits for room when the queue is full
//
// PROBLEM 3 - unbounded memory, which is what we set out to fix.
//
//   If the producer outruns the disk, blocks pile up and we are back to holding
//   the whole table.
//
//   Solution: queue_depth. submit() does not return until there is room. The
//   producer can never run more than queue_depth blocks ahead of the disk, so
//   peak memory is bounded by
//
//       queue_depth  x  scatterers  x  reflections-per-block  x  16 bytes
//
//   instead of by the size of the whole table. submit() blocking is not an
//   inefficiency; it IS the memory limit.
//
// ORDER
//
//   The file needs reflections in order, but the producer may finish blocks out
//   of order. So each block carries an id, pending_ is an ordered map keyed by
//   it, and the consumer only writes the block whose id equals next_id_. An
//   early block simply waits its turn. Ids must be 0, 1, 2, ... with no gaps,
//   or the consumer would wait forever for one that never arrives - finish()
//   reports that as an error instead of hanging.
//
// ERRORS ON THE WRITER THREAD
//
//   An exception cannot travel from the writer thread to the producer by
//   itself, and one that escapes a thread terminates the program. So drain()
//   catches everything, stores it, and the next submit() or finish() rethrows
//   it on the producer's thread where the caller can handle it.
//
// LIFETIME WARNING
//
//   The writer thread runs until finish() joins it, and everything it touches
//   must outlive that call. Blocks are therefore MOVED into the queue and owned
//   by it. Never hand over a pointer or reference to a buffer the producer
//   intends to reuse for the next block.
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

    // Opens the file, writes the header and scatterer list, then starts the
    // writer thread. From here until finish(), that thread owns the file.
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
    // full-sized block has to be built merely to write a header.
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

    // finish() is normally called explicitly, so that errors surface. This is
    // only a backstop for an early return or an exception on the producer side.
    ~tsc_stream_writer() { try { finish(); } catch (...) { } }

    tsc_stream_writer(const tsc_stream_writer&) = delete;
    tsc_stream_writer& operator=(const tsc_stream_writer&) = delete;

    // Hand one block of reflections to the writer.
    //
    // id must run 0, 1, 2, ... overall, though blocks may arrive in any order.
    // sf_block is [scatterer][reflection within the block]; idx is
    // [dimension][same]. Both are moved from: after this call the caller's
    // copies are empty, deliberately - the writer thread owns them now and the
    // producer must not touch them again.
    //
    // This call WAITS while the queue is full. That wait is the memory bound.
    void submit(const std::size_t id, index_block idx, cvec2 sf_block)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        // Sleep until there is room - or until the writer has failed, so we can
        // report that instead of waiting for space that will never appear.
        space_.wait(lock, [this] { return pending_.size() < queue_depth_ || failed_; });
        rethrow_if_failed();
        pending_.emplace(id, item{ std::move(idx), std::move(sf_block) });
        lock.unlock();       // release first, so the woken writer does not
        work_.notify_one();  // immediately block on a mutex we still hold
    }

    // Waits for every submitted block to be written, then joins the writer
    // thread and rethrows anything that failed on it. Safe to call twice.
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

    // The writer thread's entire life: take the next block in order, write it,
    // repeat. Exits once finish() has been called and nothing is left.
    void drain()
    {
        try
        {
            for (;;)
            {
                item next;
                {
                    std::unique_lock<std::mutex> lock(mutex_);
                    // Sleep until the block we need next is present, or we are
                    // told to stop. Waiting on next_id_ specifically - rather
                    // than on "the queue is non-empty" - is what keeps the file
                    // in order when blocks arrive out of order.
                    work_.wait(lock, [this] {
                        return done_ || pending_.count(next_id_) != 0;
                    });
                    const auto found = pending_.find(next_id_);
                    if (found == pending_.end())
                    {
                        // Woken, but the block we need is not here.
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

                // The write happens OUTSIDE the mutex. That is the entire point
                // of this class: the producer computes the next block while this
                // one goes to disk. Only this thread writes, so the file needs
                // no lock of its own.
                block_type::write_tscb_reflection_block(out_, next.idx, next.sf);
            }
        }
        catch (...)
        {
            // An exception must not escape a thread - that terminates the
            // program. Park it here; submit()/finish() rethrow it on the
            // caller's thread. Waking space_ releases any producer blocked
            // waiting for room that will now never be freed.
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
