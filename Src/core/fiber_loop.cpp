#include "fiber_loop.h"
#include <algorithm>
#include <vector>
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <sys/mman.h>
#include <ucontext.h>
#endif

namespace {

//Reserved per fiber, committed as touched; the field evaluators keep their scratch thread_local
constexpr size_t fiber_stack = 512 * 1024;

struct fiber;
struct scheduler {
	const std::function<void(int)> *body = nullptr;
	const std::function<void()> *round_end = nullptr;
	std::atomic<int> *next = nullptr;
	int n = 0;
	fiber *cur = nullptr;
#ifdef _WIN32
	void *self = nullptr;
#else
	ucontext_t self{};
#endif
};
struct fiber {
	scheduler *s = nullptr;
	bool done = false;
#ifdef _WIN32
	void *h = nullptr;
#else
	ucontext_t uc{};
	void *stack = nullptr;
#endif
};
thread_local scheduler *t_sched = nullptr;

void to_scheduler(fiber *f)
{
#ifdef _WIN32
	SwitchToFiber(f->s->self);
#else
	swapcontext(&f->uc, &f->s->self);
#endif
}

void fiber_main(fiber *f)
{
	for (int i; (i = f->s->next->fetch_add(1, std::memory_order_relaxed)) < f->s->n;) (*f->s->body)(i);
	f->done = true;
	to_scheduler(f); //never resumed
}

#ifdef _WIN32
void CALLBACK fiber_entry(void *p) { fiber_main(static_cast<fiber *>(p)); }
#else
// ponytail: glibc's swapcontext saves the signal mask with a syscall, ~1 us per switch; a
// hand-written switch if that ever shows next to the device time
void fiber_entry() { fiber_main(t_sched->cur); }
#endif

bool make(fiber &f)
{
#ifdef _WIN32
	f.h = CreateFiberEx(64 * 1024, fiber_stack, FIBER_FLAG_FLOAT_SWITCH, fiber_entry, &f);
	return f.h != nullptr;
#else
	f.stack = mmap(nullptr, fiber_stack, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_STACK, -1, 0);
	if (f.stack == MAP_FAILED) { f.stack = nullptr; return false; }
	getcontext(&f.uc);
	f.uc.uc_stack.ss_sp = f.stack;
	f.uc.uc_stack.ss_size = fiber_stack;
	f.uc.uc_link = nullptr;
	makecontext(&f.uc, fiber_entry, 0);
	return true;
#endif
}

void destroy(fiber &f)
{
#ifdef _WIN32
	if (f.h) DeleteFiber(f.h);
#else
	if (f.stack) munmap(f.stack, fiber_stack);
#endif
}

} //namespace

void fiber_yield()
{
	if (t_sched->cur) to_scheduler(t_sched->cur);
	else (*t_sched->round_end)();
}

void fiber_loop(const int nfib, std::atomic<int> &next, const int n, const std::function<void(int)> &body, const std::function<void()> &round_end)
{
	scheduler s;
	s.body = &body; s.round_end = &round_end; s.next = &next; s.n = n;
	scheduler *outer = t_sched;
	t_sched = &s;
#ifdef _WIN32
	const bool was_fiber = IsThreadAFiber() != FALSE;
	s.self = was_fiber ? GetCurrentFiber() : ConvertThreadToFiberEx(nullptr, FIBER_FLAG_FLOAT_SWITCH);
#endif
	//No more fibers than indices left; a fiber that cannot be made just leaves its share to the others
	std::vector<fiber> pool(static_cast<size_t>(std::max(1, std::min(nfib, n - next.load()))));
	std::vector<fiber *> live;
	for (fiber &f : pool) {
		f.s = &s;
		if (make(f)) live.push_back(&f);
	}
	if (live.empty()) //no fiber at all: the body runs here and a park is a round of one
		for (int i; (i = next.fetch_add(1, std::memory_order_relaxed)) < n;) body(i);
	while (!live.empty()) {
		for (fiber *f : live) {
			s.cur = f;
#ifdef _WIN32
			SwitchToFiber(f->h);
#else
			swapcontext(&s.self, &f->uc);
#endif
		}
		live.erase(std::remove_if(live.begin(), live.end(), [](const fiber *f) { return f->done; }), live.end());
		round_end();
	}
	for (fiber &f : pool) destroy(f);
#ifdef _WIN32
	if (!was_fiber) ConvertFiberToThread();
#endif
	t_sched = outer;
}
