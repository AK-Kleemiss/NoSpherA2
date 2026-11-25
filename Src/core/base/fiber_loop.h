#pragma once
#include <atomic>
#include <functional>

//Cooperative fibers on the calling thread, so one host thread can keep many field evaluations in
//flight for a batched device. fiber_loop runs body(i) for every i it takes from next (below n) on
//up to nfib fibers; a body parks with fiber_yield(), and once every live fiber has run to its next
//park (or finished) round_end() runs and the parked fibers resume. Fibers never leave the thread
//that started them, so thread_locals stay valid across a park.
//A body must not hold a lock or sit inside an OpenMP construct across a fiber_yield, and must not
//throw out of body.
void fiber_loop(int nfib, std::atomic<int> &next, int n, const std::function<void(int)> &body, const std::function<void()> &round_end);
void fiber_yield();
