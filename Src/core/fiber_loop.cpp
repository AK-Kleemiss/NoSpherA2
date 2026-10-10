#include "fiber_loop.h"
#include <algorithm>
#include <cstdint>
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
#ifndef MAP_STACK //Linux only; macOS stacks need no flag
#define MAP_STACK 0
#endif
#if defined(__linux__) && defined(__x86_64__)
#define NOS_FIBER_ASM 1
#else
#include <ucontext.h>
#endif
#endif

#ifdef NOS_FIBER_ASM
//System V x86-64 switch: callee-saved registers, MXCSR and the x87 control word go on the old stack, the stack pointers
//swap, and the same come off the new one. glibc's swapcontext also saves the signal mask, a syscall per switch.
// ponytail: no CET shadow stack; a binary that runs with one enabled needs swapcontext back for this file
extern "C" __attribute__((visibility("hidden"))) void nos_fiber_switch(void **save, void *load);
extern "C" __attribute__((visibility("hidden"))) void nos_fiber_start();
asm(R"(
	.text
	.p2align 4
	.globl nos_fiber_switch
	.hidden nos_fiber_switch
	.type nos_fiber_switch, @function
nos_fiber_switch:
	pushq %rbp
	pushq %rbx
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15
	subq $8, %rsp
	stmxcsr (%rsp)
	fnstcw 4(%rsp)
	movq %rsp, (%rdi)
	movq %rsi, %rsp
	ldmxcsr (%rsp)
	fldcw 4(%rsp)
	addq $8, %rsp
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbx
	popq %rbp
	ret
	.size nos_fiber_switch, .-nos_fiber_switch
	.p2align 4
	.globl nos_fiber_start
	.hidden nos_fiber_start
	.type nos_fiber_start, @function
nos_fiber_start:
	movq %r12, %rdi
	callq *%r13
	ud2
	.size nos_fiber_start, .-nos_fiber_start
)");
#endif

namespace {

//Reserved per fiber, committed as touched; the field evaluators keep their scratch thread_local
constexpr size_t fiber_stack = 512 * 1024;

//One side of a switch: a Windows fiber, a saved stack pointer, or a ucontext
#ifdef _WIN32
struct ctx { void *h = nullptr; };
#elif defined(NOS_FIBER_ASM)
struct ctx { void *sp = nullptr; };
#else
struct ctx { ucontext_t uc{}; };
#endif
void switch_to(ctx &from, ctx &to)
{
#ifdef _WIN32
	(void)from;
	SwitchToFiber(to.h);
#elif defined(NOS_FIBER_ASM)
	nos_fiber_switch(&from.sp, to.sp);
#else
	swapcontext(&from.uc, &to.uc);
#endif
}

struct fiber;
struct scheduler {
	const std::function<void(int)> *body = nullptr;
	const std::function<void()> *round_end = nullptr;
	std::atomic<int> *next = nullptr;
	int n = 0;
	fiber *cur = nullptr;
	ctx self;
};
struct fiber {
	scheduler *s = nullptr;
	bool done = false;
	ctx c;
#ifndef _WIN32
	void *stack = nullptr;
#endif
};
thread_local scheduler *t_sched = nullptr;

void to_scheduler(fiber *f) { switch_to(f->c, f->s->self); }

void fiber_main(fiber *f)
{
	for (int i; (i = f->s->next->fetch_add(1, std::memory_order_relaxed)) < f->s->n;) (*f->s->body)(i);
	f->done = true;
	to_scheduler(f); //never resumed
}

#ifdef _WIN32
void CALLBACK fiber_entry(void *p) { fiber_main(static_cast<fiber *>(p)); }
#elif defined(NOS_FIBER_ASM)
void fiber_entry(fiber *f) { fiber_main(f); }
#else
void fiber_entry() { fiber_main(t_sched->cur); }
#endif

bool make(fiber &f)
{
#ifdef _WIN32
	f.c.h = CreateFiberEx(64 * 1024, fiber_stack, FIBER_FLAG_FLOAT_SWITCH, fiber_entry, &f);
	return f.c.h != nullptr;
#else
	f.stack = mmap(nullptr, fiber_stack, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_STACK, -1, 0);
	if (f.stack == MAP_FAILED) { f.stack = nullptr; return false; }
#ifdef NOS_FIBER_ASM
	//The frame nos_fiber_switch pops: control words, r15, r14, r13 = entry, r12 = fiber, rbx, rbp, then nos_fiber_start
	//as the return address, which leaves the stack 16-byte aligned at its call
	uint64_t *sp = reinterpret_cast<uint64_t *>(static_cast<char *>(f.stack) + fiber_stack) - 10;
	uint32_t mxcsr;
	uint16_t fpucw;
	asm volatile("stmxcsr %0" : "=m"(mxcsr));
	asm volatile("fnstcw %0" : "=m"(fpucw));
	sp[0] = mxcsr | static_cast<uint64_t>(fpucw) << 32;
	sp[1] = sp[2] = sp[5] = sp[6] = 0;
	sp[3] = reinterpret_cast<uint64_t>(&fiber_entry);
	sp[4] = reinterpret_cast<uint64_t>(&f);
	sp[7] = reinterpret_cast<uint64_t>(&nos_fiber_start);
	f.c.sp = sp;
#else
	getcontext(&f.c.uc);
	f.c.uc.uc_stack.ss_sp = f.stack;
	f.c.uc.uc_stack.ss_size = fiber_stack;
	f.c.uc.uc_link = nullptr;
	makecontext(&f.c.uc, fiber_entry, 0);
#endif
	return true;
#endif
}

void destroy(fiber &f)
{
#ifdef _WIN32
	if (f.c.h) DeleteFiber(f.c.h);
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
	s.self.h = was_fiber ? GetCurrentFiber() : ConvertThreadToFiberEx(nullptr, FIBER_FLAG_FLOAT_SWITCH);
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
			switch_to(s.self, f->c);
		}
		s.cur = nullptr;
		live.erase(std::remove_if(live.begin(), live.end(), [](const fiber *f) { return f->done; }), live.end());
		round_end();
	}
	for (fiber &f : pool) destroy(f);
#ifdef _WIN32
	if (!was_fiber) ConvertFiberToThread();
#endif
	t_sched = outer;
}
