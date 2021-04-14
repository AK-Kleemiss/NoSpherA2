#define __NV_MODULE_ID _36_structure_factors_compute_61_cpp1_ii_415e9718
#define __NV_CUBIN_HANDLE_STORAGE__ static
#if !defined(__CUDA_INCLUDE_COMPILER_INTERNAL_HEADERS__)
#define __CUDA_INCLUDE_COMPILER_INTERNAL_HEADERS__
#endif
#include "crt/host_runtime.h"
#include "structure_factors.fatbin.c"
extern void __device_stub__Z9KernelWFNPdS_S_S_S_PiS0_S_S_S_ii(double *, double *, double *, double *, double *, int *, int *, double *, double *, double *, int, int);
static void __device_stub__ZN3cub11EmptyKernelIvEEvv(void);
static void __nv_cudaEntityRegisterCallback(void **);
static void __sti____cudaRegisterAll(void);
#pragma section(".CRT$XCU",read)
__declspec(allocate(".CRT$XCU"))static void (*__dummy_static_init__sti____cudaRegisterAll[])(void) = {__sti____cudaRegisterAll};
void __device_stub__Z9KernelWFNPdS_S_S_S_PiS0_S_S_S_ii(
double *__par0, 
double *__par1, 
double *__par2, 
double *__par3, 
double *__par4, 
int *__par5, 
int *__par6, 
double *__par7, 
double *__par8, 
double *__par9, 
int __par10, 
int __par11)
{
__cudaLaunchPrologue(12);
__cudaSetupArgSimple(__par0, 0Ui64);
__cudaSetupArgSimple(__par1, 8Ui64);
__cudaSetupArgSimple(__par2, 16Ui64);
__cudaSetupArgSimple(__par3, 24Ui64);
__cudaSetupArgSimple(__par4, 32Ui64);
__cudaSetupArgSimple(__par5, 40Ui64);
__cudaSetupArgSimple(__par6, 48Ui64);
__cudaSetupArgSimple(__par7, 56Ui64);
__cudaSetupArgSimple(__par8, 64Ui64);
__cudaSetupArgSimple(__par9, 72Ui64);
__cudaSetupArgSimple(__par10, 80Ui64);
__cudaSetupArgSimple(__par11, 84Ui64);
__cudaLaunch(((char *)((void ( *)(double *, double *, double *, double *, double *, int *, int *, double *, double *, double *, int, int))KernelWFN)));
}
void KernelWFN( double *__cuda_0,double *__cuda_1,double *__cuda_2,double *__cuda_3,double *__cuda_4,int *__cuda_5,int *__cuda_6,double *__cuda_7,double *__cuda_8,double *__cuda_9,int __cuda_10,int __cuda_11)
{__device_stub__Z9KernelWFNPdS_S_S_S_PiS0_S_S_S_ii( __cuda_0,__cuda_1,__cuda_2,__cuda_3,__cuda_4,__cuda_5,__cuda_6,__cuda_7,__cuda_8,__cuda_9,__cuda_10,__cuda_11);
}
#line 1 "x64/Release/structure_factors.compute_50.cudafe1.stub.c"
static void __device_stub__ZN3cub11EmptyKernelIvEEvv(void)
{
__cudaLaunchPrologue(1);
__cudaLaunch(((char *)((void ( *)(void))cub::EmptyKernel<void> )));
}namespace cub{

template<> __specialization_static void __wrapper__device_stub_EmptyKernel<void>(void){__device_stub__ZN3cub11EmptyKernelIvEEvv();}}
static void __nv_cudaEntityRegisterCallback(
void **__T4)
{
__nv_dummy_param_ref(__T4);
__nv_save_fatbinhandle_for_managed_rt(__T4);
__cudaRegisterEntry(__T4, ((void ( *)(void))cub::EmptyKernel<void> ), _ZN3cub11EmptyKernelIvEEvv, (-1));
__cudaRegisterEntry(__T4, ((void ( *)(double *, double *, double *, double *, double *, int *, int *, double *, double *, double *, int, int))KernelWFN), _Z9KernelWFNPdS_S_S_S_PiS0_S_S_S_ii, (-1));
__cudaRegisterVariable(__T4, __shadow_var(gpu_nex,::gpu_nex), 0, 4Ui64, 1, 0);
__cudaRegisterVariable(__T4, __shadow_var(gpu_nmo,::gpu_nmo), 0, 4Ui64, 1, 0);
__cudaRegisterVariable(__T4, __shadow_var(gpu_type_vector,::gpu_type_vector), 0, 672Ui64, 1, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust6system6detail10sequential3seqE,::thrust::system::detail::sequential::seq), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust6system3cpp3parE,::thrust::system::cpp::par), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust8cuda_cub3parE,::thrust::cuda_cub::par), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders2_1E,::thrust::placeholders::_1), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders2_2E,::thrust::placeholders::_2), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders2_3E,::thrust::placeholders::_3), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders2_4E,::thrust::placeholders::_4), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders2_5E,::thrust::placeholders::_5), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders2_6E,::thrust::placeholders::_6), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders2_7E,::thrust::placeholders::_7), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders2_8E,::thrust::placeholders::_8), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders2_9E,::thrust::placeholders::_9), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust12placeholders3_10E,::thrust::placeholders::_10), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust3seqE,::thrust::seq), 0, 1Ui64, 0, 0);
__cudaRegisterVariable(__T4, __shadow_var(__nv_static_49__36_structure_factors_compute_61_cpp1_ii_415e9718__ZN58_INTERNAL_36_structure_factors_compute_61_cpp1_ii_415e97186thrust6deviceE,::thrust::device), 0, 1Ui64, 0, 0);
}
static void __sti____cudaRegisterAll(void)
{
____cudaRegisterLinkedBinary(__nv_cudaEntityRegisterCallback);
}
