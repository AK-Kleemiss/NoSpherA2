#include "itensor_gpu.h"
#include "throughput.h"
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>
#include <cstdio>
#include <cstdint>
#include <map>
#include <vector>
#include <algorithm>

//The I tensor's device path on Apple silicon, behind the same six entry points as the
//CUDA one so XCW.cpp does not know which it got. The GEMM is Metal Performance Shaders'
//single-precision one; phase, weighting and accumulation are the three kernels below,
//compiled from source when the device is opened so the build needs no Metal compiler.
//
//Two things differ from the CUDA path. Memory is unified, so the AO values are converted
//once into a buffer the GPU reads in place and the results are read back without a copy.
//And Metal has no double: the phase argument and the packed accumulation, which the CUDA
//path keeps in double, are done in float-float here - a value is a pair (hi, lo) carrying
//about 48 bits of mantissa - and added into the caller's double on the host. Asking for
//FP64 therefore keeps the CPU path.
//
//Every buffer offset handed to Metal is a multiple of 64 bytes: rows are padded to a
//multiple of 16 floats, which is also what lets one MPS matrix span the real and the
//imaginary halves of the weighted copy.

namespace {

constexpr int row_pad = 16;
constexpr size_t align_bytes = 64;

inline int padded(const int n) { return (n + row_pad - 1) / row_pad * row_pad; }

//The kernels. Error-free transforms need the compiler to leave the arithmetic alone, so the
//library is compiled with fast math off (see open_device).
const char* kernel_source = R"MSL(
#include <metal_stdlib>
using namespace metal;

inline float2 two_sum(const float a, const float b)
{
	const float s = a + b;
	const float bb = s - a;
	return float2(s, (a - (s - bb)) + (b - bb));
}
inline float2 two_prod(const float a, const float b)
{
	const float p = a * b;
	return float2(p, fma(a, b, -p));
}
inline float2 df_add(const float2 x, const float2 y)
{
	float2 s = two_sum(x.x, y.x);
	s.y += x.y + y.y;
	const float h = s.x + s.y;
	return float2(h, s.y - (h - s.x));
}
inline float2 df_mul(const float2 x, const float2 y)
{
	float2 p = two_prod(x.x, y.x);
	p.y += x.x * y.y + x.y * y.x;
	const float h = p.x + p.y;
	return float2(h, p.y - (h - p.x));
}

struct phase_args { float2 kx, ky, kz; uint n; };

//k arrives divided by 2pi, so t is in turns and the reduction is a rint. Coordinates and
//k are float-float pairs: k.d reaches a few hundred turns, and a plain float there would
//cost 1e-5 in the angle where the CUDA path, in double, costs nothing.
kernel void phase_kernel(device const float2* d1 [[buffer(0)]],
	device const float2* d2 [[buffer(1)]], device const float2* d3 [[buffer(2)]],
	device const float* w [[buffer(3)]], device float* pre [[buffer(4)]],
	device float* pim [[buffer(5)]], constant phase_args& a [[buffer(6)]],
	uint p [[thread_position_in_grid]])
{
	if (p >= a.n) return;
	const float2 t = df_add(df_add(df_mul(a.kx, d1[p]), df_mul(a.ky, d2[p])), df_mul(a.kz, d3[p]));
	const float frac = (t.x - rint(t.x)) + t.y;
	const float ang = 2.0f * M_PI_F * frac;
	pre[p] = w[p] * cos(ang);
	pim[p] = w[p] * sin(ang);
}

struct weight_args { uint na, np, ld, point_base; ulong ao_off; };

//ao is row-major n_active x ld per block; the two weighted copies are written with the same
//leading dimension, the imaginary one straight after the real one.
kernel void weight_kernel(device const float* ao [[buffer(0)]],
	device const float* pre [[buffer(1)]], device const float* pim [[buffer(2)]],
	device float* wre [[buffer(3)]], device float* wim [[buffer(4)]],
	constant weight_args& a [[buffer(5)]], uint2 id [[thread_position_in_grid]])
{
	const uint p = id.x, i = id.y;
	if (p >= a.np || i >= a.na) return;
	const ulong e = (ulong)i * a.ld + p;
	const float v = ao[a.ao_off + e];
	wre[e] = v * pre[a.point_base + p];
	wim[e] = v * pim[a.point_base + p];
}

struct acc_args { uint na, nmo, ldc, aos_off; float fre, fim; };

//C is row-major 2 n_active x n_active: the real rows first, the imaginary ones below. It is
//symmetric, so only i <= j is read. Distinct (i,j) map to distinct stored indices and the
//dispatches of one command buffer run in order, so the read-modify-write needs no atomics.
kernel void accumulate_kernel(device const int* aos [[buffer(0)]],
	device const int* compact [[buffer(1)]], device const float* c [[buffer(2)]],
	device float* I_re_hi [[buffer(3)]], device float* I_re_lo [[buffer(4)]],
	device float* I_im_hi [[buffer(5)]], device float* I_im_lo [[buffer(6)]],
	constant acc_args& a [[buffer(7)]], uint2 id [[thread_position_in_grid]])
{
	const uint j = id.x, i = id.y;
	if (i >= a.na || j >= a.na || j < i) return;
	const long mu = aos[a.aos_off + i];
	const long nu = aos[a.aos_off + j];
	const int t = compact[mu * (long)a.nmo + nu];
	if (t < 0) return;
	const float re = c[(ulong)j * a.ldc + i];
	const float im = c[(ulong)(a.na + j) * a.ldc + i];
	const float2 dre = df_add(two_prod(re, a.fre), two_prod(-im, a.fim));
	const float2 dim = df_add(two_prod(re, a.fim), two_prod(im, a.fre));
	const float2 sre = df_add(float2(I_re_hi[t], I_re_lo[t]), dre);
	const float2 sim = df_add(float2(I_im_hi[t], I_im_lo[t]), dim);
	I_re_hi[t] = sre.x; I_re_lo[t] = sre.y;
	I_im_hi[t] = sim.x; I_im_lo[t] = sim.y;
}
)MSL";

//Host mirrors of the argument structs. float2 aligns to 8 in the shader, which these
//layouts reproduce by hand.
struct phase_args { float kx[2], ky[2], kz[2]; uint32_t n, pad; };
struct weight_args { uint32_t na, np, ld, point_base; uint64_t ao_off; };
struct acc_args { uint32_t na, nmo, ldc, aos_off; float fre, fim; };

inline void split(const double v, float* hi, float* lo)
{
	*hi = (float)v;
	*lo = (float)(v - (double)*hi);
}

struct Device {
	id<MTLDevice> dev = nil;
	id<MTLCommandQueue> queue = nil;
	id<MTLComputePipelineState> phase = nil, weight = nil, accumulate = nil;
	bool opened = false;
	bool ok = false;
};

Device& device()
{
	static Device d;
	return d;
}

bool open_device()
{
	Device& d = device();
	if (d.opened) return d.ok;
	d.opened = true;
	@autoreleasepool {
		d.dev = MTLCreateSystemDefaultDevice();
		if (d.dev == nil) return false;
		//Shared buffers are the whole design; a discrete card behind PCIe would read the
		//AO values over the bus on every reflection.
		if (!d.dev.hasUnifiedMemory) return false;
		MTLCompileOptions* opts = [MTLCompileOptions new];
		if (@available(macOS 15.0, *)) {
			opts.mathMode = MTLMathModeSafe;
			opts.mathFloatingPointFunctions = MTLMathFloatingPointFunctionsPrecise;
		} else {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
			opts.fastMathEnabled = NO;
#pragma clang diagnostic pop
		}
		NSError* err = nil;
		id<MTLLibrary> lib = [d.dev newLibraryWithSource:[NSString stringWithUTF8String:kernel_source]
			options:opts error:&err];
		if (lib == nil) {
			std::fprintf(stderr, "NoSpherA2 I tensor Metal: shader compilation failed: %s\n",
				err ? err.localizedDescription.UTF8String : "unknown error");
			return false;
		}
		auto pipeline = [&](const char* name) -> id<MTLComputePipelineState> {
			id<MTLFunction> f = [lib newFunctionWithName:[NSString stringWithUTF8String:name]];
			if (f == nil) return nil;
			NSError* e = nil;
			id<MTLComputePipelineState> p = [d.dev newComputePipelineStateWithFunction:f error:&e];
			if (p == nil)
				std::fprintf(stderr, "NoSpherA2 I tensor Metal: %s: %s\n", name,
					e ? e.localizedDescription.UTF8String : "no pipeline");
			return p;
		};
		d.phase = pipeline("phase_kernel");
		d.weight = pipeline("weight_kernel");
		d.accumulate = pipeline("accumulate_kernel");
		d.queue = [d.dev newCommandQueue];
		d.ok = d.phase != nil && d.weight != nil && d.accumulate != nil && d.queue != nil;
	}
	return d.ok;
}

struct Block {
	int na = 0, np = 0, ld = 0, ldc = 0, grid = 0, point_start = 0;
	uint64_t ao_off = 0;   //floats, into the padded AO buffer
	uint32_t aos_off = 0;
	MPSMatrix* A = nil;    //n_active x np over the AO buffer
	MPSMatrix* W = nil;    //2 n_active x np over the weighted scratch
	MPSMatrix* C = nil;    //2 n_active x n_active over the result scratch
	MPSMatrixMultiplication* gemm = nil;
};

struct State {
	bool ready = false;
	int nmo = 0, packed = 0, n_grids = 0;
	double issued_flops = 0.0;
	long long n_points = 0;
	id<MTLBuffer> ao = nil, aos = nil, compact = nil;
	id<MTLBuffer> d1 = nil, d2 = nil, d3 = nil, w = nil;
	id<MTLBuffer> phase_re = nil, phase_im = nil;
	id<MTLBuffer> wri = nil, cri = nil;
	//Per slot: four packed float arrays, re hi, re lo, im hi, im lo, back to back
	id<MTLBuffer> I[2] = { nil, nil };
	id<MTLCommandBuffer> pending[2] = { nil, nil };
	std::vector<Block> blocks;
	std::vector<int> grid_off;
	//One MPS kernel per distinct (n_active, points) shape
	std::map<std::pair<int, int>, MPSMatrixMultiplication*> gemms;
};

State& state()
{
	static State s;
	return s;
}

id<MTLBuffer> make_buffer(const size_t bytes)
{
	return [device().dev newBufferWithLength:std::max<size_t>(bytes, align_bytes)
		options:MTLResourceStorageModeShared];
}

MPSMatrix* matrix(id<MTLBuffer> buf, const size_t offset_bytes, const int rows, const int cols, const int ld)
{
	MPSMatrixDescriptor* desc = [MPSMatrixDescriptor matrixDescriptorWithRows:rows columns:cols
		rowBytes:(size_t)ld * sizeof(float) dataType:MPSDataTypeFloat32];
	return [[MPSMatrix alloc] initWithBuffer:buf offset:offset_bytes descriptor:desc];
}

bool init_impl(const itensor_gpu_layout& L)
{
	State& s = state();
	int max_na = 0, max_ld = 0, max_ldc = 0;
	uint64_t ao_len = 0;
	for (int b = 0; b < L.n_blocks; b++) {
		const int na = L.blk_n_active[b];
		const int ld = padded(L.blk_point_count[b]);
		max_na = std::max(max_na, na);
		max_ld = std::max(max_ld, ld);
		max_ldc = std::max(max_ldc, padded(na));
		ao_len += (uint64_t)na * ld;
		s.issued_flops += throughput::flops_gemm(na, 2.0 * na, L.blk_point_count[b]);
	}
	if (throughput::enabled()) {
		int min_na = L.n_blocks ? L.blk_n_active[0] : 0;
		int min_np = L.n_blocks ? L.blk_point_count[0] : 0, max_np = 0;
		for (int b = 0; b < L.n_blocks; b++) {
			min_na = std::min(min_na, L.blk_n_active[b]);
			min_np = std::min(min_np, L.blk_point_count[b]);
			max_np = std::max(max_np, L.blk_point_count[b]);
		}
		std::fprintf(stderr, "I tensor GPU: %d blocks, n_active %d-%d, points %d-%d\n",
			L.n_blocks, min_na, max_na, min_np, max_np);
	}
	const size_t need =
		sizeof(float) * (size_t)ao_len +
		sizeof(int) * (size_t)L.aos_all_len +
		sizeof(int) * (size_t)L.nmo * L.nmo +
		sizeof(float) * 7 * (size_t)L.n_points +
		sizeof(float) * 2 * (size_t)max_na * max_ld +
		sizeof(float) * 2 * (size_t)max_na * max_ldc +
		sizeof(float) * 8 * (size_t)L.packed;
	const uint64_t limit = device().dev.recommendedMaxWorkingSetSize;
	if (need + (1u << 28) > limit) {
		std::fprintf(stderr, "I tensor Metal: %.1f GB needed, %.1f GB available to the GPU\n",
			need / 1.0e9, limit / 1.0e9);
		return false;
	}

	s.ao = make_buffer(sizeof(float) * ao_len);
	s.aos = make_buffer(sizeof(int) * L.aos_all_len);
	s.compact = make_buffer(sizeof(int) * (size_t)L.nmo * L.nmo);
	s.d1 = make_buffer(sizeof(float) * 2 * L.n_points);
	s.d2 = make_buffer(sizeof(float) * 2 * L.n_points);
	s.d3 = make_buffer(sizeof(float) * 2 * L.n_points);
	s.w = make_buffer(sizeof(float) * L.n_points);
	s.phase_re = make_buffer(sizeof(float) * L.n_points);
	s.phase_im = make_buffer(sizeof(float) * L.n_points);
	s.wri = make_buffer(sizeof(float) * 2 * (size_t)max_na * max_ld);
	s.cri = make_buffer(sizeof(float) * 2 * (size_t)max_na * max_ldc);
	for (int i = 0; i < 2; i++) s.I[i] = make_buffer(sizeof(float) * 4 * (size_t)L.packed);
	for (id<MTLBuffer> b : { s.ao, s.aos, s.compact, s.d1, s.d2, s.d3, s.w, s.phase_re, s.phase_im,
		s.wri, s.cri, s.I[0], s.I[1] })
		if (b == nil) return false;

	//The AO values, narrowed and re-laid with padded rows, straight into the shared buffer
	float* ao = (float*)s.ao.contents;
	s.blocks.resize(L.n_blocks);
	uint64_t off = 0;
	for (int b = 0; b < L.n_blocks; b++) {
		Block& blk = s.blocks[b];
		blk.na = L.blk_n_active[b];
		blk.np = L.blk_point_count[b];
		blk.ld = padded(blk.np);
		blk.ldc = padded(blk.na);
		blk.grid = L.blk_grid[b];
		blk.point_start = L.blk_point_start[b];
		blk.ao_off = off;
		blk.aos_off = (uint32_t)L.blk_aos_off[b];
		const double* src = L.ao_all + L.blk_ao_off[b];
		for (int i = 0; i < blk.na; i++) {
			float* row = ao + off + (uint64_t)i * blk.ld;
			for (int p = 0; p < blk.np; p++) row[p] = (float)src[(uint64_t)i * blk.np + p];
			for (int p = blk.np; p < blk.ld; p++) row[p] = 0.0f;
		}
		off += (uint64_t)blk.na * blk.ld;
		blk.A = matrix(s.ao, sizeof(float) * blk.ao_off, blk.na, blk.np, blk.ld);
		blk.W = matrix(s.wri, 0, 2 * blk.na, blk.np, blk.ld);
		blk.C = matrix(s.cri, 0, 2 * blk.na, blk.na, blk.ldc);
		const std::pair<int, int> shape(blk.na, blk.np);
		auto it = s.gemms.find(shape);
		if (it == s.gemms.end()) {
			//C = W A^T, row-major throughout: (2na x np) times (np x na)
			MPSMatrixMultiplication* g = [[MPSMatrixMultiplication alloc] initWithDevice:device().dev
				transposeLeft:NO transposeRight:YES resultRows:2 * blk.na resultColumns:blk.na
				interiorColumns:blk.np alpha:1.0 beta:0.0];
			it = s.gemms.emplace(shape, g).first;
		}
		blk.gemm = it->second;
	}
	std::copy(L.aos_all, L.aos_all + L.aos_all_len, (int*)s.aos.contents);
	std::copy(L.compact, L.compact + (size_t)L.nmo * L.nmo, (int*)s.compact.contents);
	float* d1 = (float*)s.d1.contents;
	float* d2 = (float*)s.d2.contents;
	float* d3 = (float*)s.d3.contents;
	float* w = (float*)s.w.contents;
	for (long long p = 0; p < L.n_points; p++) {
		split(L.d1[p], d1 + 2 * p, d1 + 2 * p + 1);
		split(L.d2[p], d2 + 2 * p, d2 + 2 * p + 1);
		split(L.d3[p], d3 + 2 * p, d3 + 2 * p + 1);
		w[p] = (float)L.weights[p];
	}
	s.grid_off.assign(L.grid_point_off, L.grid_point_off + L.n_grids + 1);
	s.nmo = L.nmo; s.packed = L.packed; s.n_grids = L.n_grids; s.n_points = L.n_points;
	s.ready = true;
	return true;
}

bool submit_impl(const int slot, const int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors)
{
	State& s = state();
	Device& d = device();
	if (!s.ready || slot < 0 || slot > 1) return false;
	@autoreleasepool {
		id<MTLCommandBuffer> cb = [d.queue commandBuffer];
		if (cb == nil) return false;
		{
			id<MTLBlitCommandEncoder> blit = [cb blitCommandEncoder];
			[blit fillBuffer:s.I[slot] range:NSMakeRange(0, sizeof(float) * 4 * (size_t)s.packed) value:0];
			[blit endEncoding];
		}
		const size_t I_bytes = sizeof(float) * (size_t)s.packed;
		for (int sy = 0; sy < num_syms; sy++) {
			phase_args pa;
			split(kx[sy] * 0.15915494309189533576888376337251, &pa.kx[0], &pa.kx[1]);
			split(ky[sy] * 0.15915494309189533576888376337251, &pa.ky[0], &pa.ky[1]);
			split(kz[sy] * 0.15915494309189533576888376337251, &pa.kz[0], &pa.kz[1]);
			pa.n = (uint32_t)s.n_points;
			pa.pad = 0;
			{
				id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
				[enc setComputePipelineState:d.phase];
				[enc setBuffer:s.d1 offset:0 atIndex:0];
				[enc setBuffer:s.d2 offset:0 atIndex:1];
				[enc setBuffer:s.d3 offset:0 atIndex:2];
				[enc setBuffer:s.w offset:0 atIndex:3];
				[enc setBuffer:s.phase_re offset:0 atIndex:4];
				[enc setBuffer:s.phase_im offset:0 atIndex:5];
				[enc setBytes:&pa length:sizeof(pa) atIndex:6];
				[enc dispatchThreads:MTLSizeMake((NSUInteger)s.n_points, 1, 1)
					threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
				[enc endEncoding];
			}
			for (const Block& blk : s.blocks) {
				const std::complex<double> f = factors[(size_t)sy * s.n_grids + blk.grid];
				if (f.real() == 0.0 && f.imag() == 0.0) continue;
				const size_t half = sizeof(float) * (size_t)blk.na * blk.ld;
				weight_args wa{ (uint32_t)blk.na, (uint32_t)blk.np, (uint32_t)blk.ld,
					(uint32_t)(s.grid_off[blk.grid] + blk.point_start), blk.ao_off };
				{
					id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
					[enc setComputePipelineState:d.weight];
					[enc setBuffer:s.ao offset:0 atIndex:0];
					[enc setBuffer:s.phase_re offset:0 atIndex:1];
					[enc setBuffer:s.phase_im offset:0 atIndex:2];
					[enc setBuffer:s.wri offset:0 atIndex:3];
					[enc setBuffer:s.wri offset:half atIndex:4];
					[enc setBytes:&wa length:sizeof(wa) atIndex:5];
					[enc dispatchThreads:MTLSizeMake(blk.np, blk.na, 1)
						threadsPerThreadgroup:MTLSizeMake(32, 8, 1)];
					[enc endEncoding];
				}
				[blk.gemm encodeToCommandBuffer:cb leftMatrix:blk.W rightMatrix:blk.A resultMatrix:blk.C];
				acc_args aa{ (uint32_t)blk.na, (uint32_t)s.nmo, (uint32_t)blk.ldc, blk.aos_off,
					(float)f.real(), (float)f.imag() };
				{
					id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
					[enc setComputePipelineState:d.accumulate];
					[enc setBuffer:s.aos offset:0 atIndex:0];
					[enc setBuffer:s.compact offset:0 atIndex:1];
					[enc setBuffer:s.cri offset:0 atIndex:2];
					for (int q = 0; q < 4; q++) [enc setBuffer:s.I[slot] offset:q * I_bytes atIndex:3 + q];
					[enc setBytes:&aa length:sizeof(aa) atIndex:7];
					[enc dispatchThreads:MTLSizeMake(blk.na, blk.na, 1)
						threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
					[enc endEncoding];
				}
			}
		}
		[cb commit];
		s.pending[slot] = cb;
	}
	return true;
}

bool collect_impl(const int slot, std::complex<double>* I_r)
{
	State& s = state();
	if (!s.ready || slot < 0 || slot > 1 || s.pending[slot] == nil) return false;
	@autoreleasepool {
		id<MTLCommandBuffer> cb = s.pending[slot];
		[cb waitUntilCompleted];
		s.pending[slot] = nil;
		if (cb.status != MTLCommandBufferStatusCompleted) {
			std::fprintf(stderr, "NoSpherA2 I tensor Metal: command buffer failed: %s\n",
				cb.error ? cb.error.localizedDescription.UTF8String : "unknown error");
			return false;
		}
	}
	const float* re_hi = (const float*)s.I[slot].contents;
	const float* re_lo = re_hi + s.packed;
	const float* im_hi = re_lo + s.packed;
	const float* im_lo = im_hi + s.packed;
	for (int i = 0; i < s.packed; i++)
		I_r[i] += std::complex<double>((double)re_hi[i] + (double)re_lo[i],
			(double)im_hi[i] + (double)im_lo[i]);
	return true;
}

} //namespace

bool itensor_gpu_available() { return open_device(); }

const char* itensor_gpu_gemm_name() { return "Metal Performance Shaders"; }

double itensor_gpu_issued_flops() { return state().issued_flops; }

bool itensor_gpu_init(const itensor_gpu_layout& L, const sf_precision prec, const bool)
{
	itensor_gpu_free();
	if (!itensor_gpu_available()) return false;
	if (prec == sf_precision::FP64) {
		std::fprintf(stderr, "I tensor Metal: no double precision on this device, keeping the CPU path\n");
		return false;
	}
	bool ok = false;
	@autoreleasepool { ok = init_impl(L); }
	if (!ok) itensor_gpu_free();
	return ok;
}

//The Metal engine contracts one reflection at a time, so a batch is one reflection here
int itensor_gpu_batch(const int) { return 1; }

bool itensor_gpu_submit(const int slot, const int n_refl, const int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors)
{
	if (n_refl != 1) return false;
	return submit_impl(slot, num_syms, kx, ky, kz, factors);
}

bool itensor_gpu_collect(const int slot, std::complex<double>* I_r, const long long)
{
	return collect_impl(slot, I_r);
}

void itensor_gpu_free()
{
	State& s = state();
	@autoreleasepool {
		for (int i = 0; i < 2; i++)
			if (s.pending[i] != nil) [s.pending[i] waitUntilCompleted];
		s = State{};
	}
}
