#define __NV_MODULE_ID _25_structure_factors_cpp1_ii_415e9718
#define __NV_CUBIN_HANDLE_STORAGE__ static
#if !defined(__CUDA_INCLUDE_COMPILER_INTERNAL_HEADERS__)
#define __CUDA_INCLUDE_COMPILER_INTERNAL_HEADERS__
#endif
#include "crt/host_runtime.h"
#include "structure_factors.fatbin.c"
extern void __device_stub__Z13gpu_make_gridiPKdS0_S0_S0_S0_S0_S0_PKiiPdS3_S3_S3_S3_i(const int, const double *, const double *, const double *, const double *, const double *, const double *, const double *, const int *, const int, double *, double *, double *, double *, double *, const int);
extern void __device_stub__Z40gpu_linear_interpolate_spherical_densityPKdS0_ibiS0_S0_S0_PdS1_S0_S0_S0_iddi(const double *, const double *, const int, const bool, const int, const double *, const double *, const double *, double *, double *, const double *, const double *, const double *, const int, const double, const double, const int);
extern void __device_stub__Z17gpu_apply_weightsiPKdS0_PdS0_i(const int, const double *, const double *, double *, const double *, const int);
extern void __device_stub__Z18gpu_count_non_zeroiPKdPi(const int, const double *, int *);
extern void __device_stub__Z16gpu_calc_chargesPKdS0_S0_S0_S0_iiPd(const double *, const double *, const double *, const double *, const double *, const int, const int, double *);
extern void __device_stub__Z13gpu_calc_densPdPKdS1_S1_S1_S1_S1_PKiS3_S1_S1_S1_(double *, const double *, const double *, const double *, const double *, const double *, const double *, const int *, const int *, const double *, const double *, const double *);
static void __nv_cudaEntityRegisterCallback(void **);
static void __sti____cudaRegisterAll(void);
#pragma section(".CRT$XCU",read)
__declspec(allocate(".CRT$XCU"))static void (*__dummy_static_init__sti____cudaRegisterAll[])(void) = {__sti____cudaRegisterAll};
void __device_stub__Z13gpu_make_gridiPKdS0_S0_S0_S0_S0_S0_PKiiPdS3_S3_S3_S3_i(
const int __par0, 
const double *__par1, 
const double *__par2, 
const double *__par3, 
const double *__par4, 
const double *__par5, 
const double *__par6, 
const double *__par7, 
const int *__par8, 
const int __par9, 
double *__par10, 
double *__par11, 
double *__par12, 
double *__par13, 
double *__par14, 
const int __par15)
{
__cudaLaunchPrologue(16);
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
__cudaSetupArgSimple(__par11, 88Ui64);
__cudaSetupArgSimple(__par12, 96Ui64);
__cudaSetupArgSimple(__par13, 104Ui64);
__cudaSetupArgSimple(__par14, 112Ui64);
__cudaSetupArgSimple(__par15, 120Ui64);
__cudaLaunch(((char *)((void ( *)(const int, const double *, const double *, const double *, const double *, const double *, const double *, const double *, const int *, const int, double *, double *, double *, double *, double *, const int))gpu_make_grid)));
}
void gpu_make_grid( const int __cuda_0,const double *__cuda_1,const double *__cuda_2,const double *__cuda_3,const double *__cuda_4,const double *__cuda_5,const double *__cuda_6,const double *__cuda_7,const int *__cuda_8,const int __cuda_9,double *__cuda_10,double *__cuda_11,double *__cuda_12,double *__cuda_13,double *__cuda_14,const int __cuda_15)
{__device_stub__Z13gpu_make_gridiPKdS0_S0_S0_S0_S0_S0_PKiiPdS3_S3_S3_S3_i( __cuda_0,__cuda_1,__cuda_2,__cuda_3,__cuda_4,__cuda_5,__cuda_6,__cuda_7,__cuda_8,__cuda_9,__cuda_10,__cuda_11,__cuda_12,__cuda_13,__cuda_14,__cuda_15);
}
#line 1 "x64/Release/structure_factors.cudafe1.stub.c"
void __device_stub__Z40gpu_linear_interpolate_spherical_densityPKdS0_ibiS0_S0_S0_PdS1_S0_S0_S0_iddi(
const double *__par0, 
const double *__par1, 
const int __par2, 
const bool __par3, 
const int __par4, 
const double *__par5, 
const double *__par6, 
const double *__par7, 
double *__par8, 
double *__par9, 
const double *__par10, 
const double *__par11, 
const double *__par12, 
const int __par13, 
const double __par14, 
const double __par15, 
const int __par16)
{
__cudaLaunchPrologue(17);
__cudaSetupArgSimple(__par0, 0Ui64);
__cudaSetupArgSimple(__par1, 8Ui64);
__cudaSetupArgSimple(__par2, 16Ui64);
__cudaSetupArgSimple(__par3, 20Ui64);
__cudaSetupArgSimple(__par4, 24Ui64);
__cudaSetupArgSimple(__par5, 32Ui64);
__cudaSetupArgSimple(__par6, 40Ui64);
__cudaSetupArgSimple(__par7, 48Ui64);
__cudaSetupArgSimple(__par8, 56Ui64);
__cudaSetupArgSimple(__par9, 64Ui64);
__cudaSetupArgSimple(__par10, 72Ui64);
__cudaSetupArgSimple(__par11, 80Ui64);
__cudaSetupArgSimple(__par12, 88Ui64);
__cudaSetupArgSimple(__par13, 96Ui64);
__cudaSetupArgSimple(__par14, 104Ui64);
__cudaSetupArgSimple(__par15, 112Ui64);
__cudaSetupArgSimple(__par16, 120Ui64);
__cudaLaunch(((char *)((void ( *)(const double *, const double *, const int, const bool, const int, const double *, const double *, const double *, double *, double *, const double *, const double *, const double *, const int, const double, const double, const int))gpu_linear_interpolate_spherical_density)));
}
void gpu_linear_interpolate_spherical_density( const double *__cuda_0,const double *__cuda_1,const int __cuda_2,const bool __cuda_3,const int __cuda_4,const double *__cuda_5,const double *__cuda_6,const double *__cuda_7,double *__cuda_8,double *__cuda_9,const double *__cuda_10,const double *__cuda_11,const double *__cuda_12,const int __cuda_13,const double __cuda_14,const double __cuda_15,const int __cuda_16)
{__device_stub__Z40gpu_linear_interpolate_spherical_densityPKdS0_ibiS0_S0_S0_PdS1_S0_S0_S0_iddi( __cuda_0,__cuda_1,__cuda_2,__cuda_3,__cuda_4,__cuda_5,__cuda_6,__cuda_7,__cuda_8,__cuda_9,__cuda_10,__cuda_11,__cuda_12,__cuda_13,__cuda_14,__cuda_15,__cuda_16);
}
#line 1 "x64/Release/structure_factors.cudafe1.stub.c"
void __device_stub__Z17gpu_apply_weightsiPKdS0_PdS0_i(
const int __par0, 
const double *__par1, 
const double *__par2, 
double *__par3, 
const double *__par4, 
const int __par5)
{
__cudaLaunchPrologue(6);
__cudaSetupArgSimple(__par0, 0Ui64);
__cudaSetupArgSimple(__par1, 8Ui64);
__cudaSetupArgSimple(__par2, 16Ui64);
__cudaSetupArgSimple(__par3, 24Ui64);
__cudaSetupArgSimple(__par4, 32Ui64);
__cudaSetupArgSimple(__par5, 40Ui64);
__cudaLaunch(((char *)((void ( *)(const int, const double *, const double *, double *, const double *, const int))gpu_apply_weights)));
}
void gpu_apply_weights( const int __cuda_0,const double *__cuda_1,const double *__cuda_2,double *__cuda_3,const double *__cuda_4,const int __cuda_5)
{__device_stub__Z17gpu_apply_weightsiPKdS0_PdS0_i( __cuda_0,__cuda_1,__cuda_2,__cuda_3,__cuda_4,__cuda_5);
}
#line 1 "x64/Release/structure_factors.cudafe1.stub.c"
void __device_stub__Z18gpu_count_non_zeroiPKdPi(
const int __par0, 
const double *__par1, 
int *__par2)
{
__cudaLaunchPrologue(3);
__cudaSetupArgSimple(__par0, 0Ui64);
__cudaSetupArgSimple(__par1, 8Ui64);
__cudaSetupArgSimple(__par2, 16Ui64);
__cudaLaunch(((char *)((void ( *)(const int, const double *, int *))gpu_count_non_zero)));
}
void gpu_count_non_zero( const int __cuda_0,const double *__cuda_1,int *__cuda_2)
{__device_stub__Z18gpu_count_non_zeroiPKdPi( __cuda_0,__cuda_1,__cuda_2);
}
#line 1 "x64/Release/structure_factors.cudafe1.stub.c"
void __device_stub__Z16gpu_calc_chargesPKdS0_S0_S0_S0_iiPd(
const double *__par0, 
const double *__par1, 
const double *__par2, 
const double *__par3, 
const double *__par4, 
const int __par5, 
const int __par6, 
double *__par7)
{
__cudaLaunchPrologue(8);
__cudaSetupArgSimple(__par0, 0Ui64);
__cudaSetupArgSimple(__par1, 8Ui64);
__cudaSetupArgSimple(__par2, 16Ui64);
__cudaSetupArgSimple(__par3, 24Ui64);
__cudaSetupArgSimple(__par4, 32Ui64);
__cudaSetupArgSimple(__par5, 40Ui64);
__cudaSetupArgSimple(__par6, 44Ui64);
__cudaSetupArgSimple(__par7, 48Ui64);
__cudaLaunch(((char *)((void ( *)(const double *, const double *, const double *, const double *, const double *, const int, const int, double *))gpu_calc_charges)));
}
void gpu_calc_charges( const double *__cuda_0,const double *__cuda_1,const double *__cuda_2,const double *__cuda_3,const double *__cuda_4,const int __cuda_5,const int __cuda_6,double *__cuda_7)
{__device_stub__Z16gpu_calc_chargesPKdS0_S0_S0_S0_iiPd( __cuda_0,__cuda_1,__cuda_2,__cuda_3,__cuda_4,__cuda_5,__cuda_6,__cuda_7);
}
#line 1 "x64/Release/structure_factors.cudafe1.stub.c"
void __device_stub__Z13gpu_calc_densPdPKdS1_S1_S1_S1_S1_PKiS3_S1_S1_S1_(
double *__par0, 
const double *__par1, 
const double *__par2, 
const double *__par3, 
const double *__par4, 
const double *__par5, 
const double *__par6, 
const int *__par7, 
const int *__par8, 
const double *__par9, 
const double *__par10, 
const double *__par11)
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
__cudaSetupArgSimple(__par11, 88Ui64);
__cudaLaunch(((char *)((void ( *)(double *, const double *, const double *, const double *, const double *, const double *, const double *, const int *, const int *, const double *, const double *, const double *))gpu_calc_dens)));
}
void gpu_calc_dens( double *__cuda_0,const double *__cuda_1,const double *__cuda_2,const double *__cuda_3,const double *__cuda_4,const double *__cuda_5,const double *__cuda_6,const int *__cuda_7,const int *__cuda_8,const double *__cuda_9,const double *__cuda_10,const double *__cuda_11)
{__device_stub__Z13gpu_calc_densPdPKdS1_S1_S1_S1_S1_PKiS3_S1_S1_S1_( __cuda_0,__cuda_1,__cuda_2,__cuda_3,__cuda_4,__cuda_5,__cuda_6,__cuda_7,__cuda_8,__cuda_9,__cuda_10,__cuda_11);
}
#line 1 "x64/Release/structure_factors.cudafe1.stub.c"
static void __nv_cudaEntityRegisterCallback(
void **__T14)
{
__nv_dummy_param_ref(__T14);
__nv_save_fatbinhandle_for_managed_rt(__T14);
__cudaRegisterEntry(__T14, ((void ( *)(double *, const double *, const double *, const double *, const double *, const double *, const double *, const int *, const int *, const double *, const double *, const double *))gpu_calc_dens), _Z13gpu_calc_densPdPKdS1_S1_S1_S1_S1_PKiS3_S1_S1_S1_, (-1));
__cudaRegisterEntry(__T14, ((void ( *)(const double *, const double *, const double *, const double *, const double *, const int, const int, double *))gpu_calc_charges), _Z16gpu_calc_chargesPKdS0_S0_S0_S0_iiPd, (-1));
__cudaRegisterEntry(__T14, ((void ( *)(const int, const double *, int *))gpu_count_non_zero), _Z18gpu_count_non_zeroiPKdPi, (-1));
__cudaRegisterEntry(__T14, ((void ( *)(const int, const double *, const double *, double *, const double *, const int))gpu_apply_weights), _Z17gpu_apply_weightsiPKdS0_PdS0_i, (-1));
__cudaRegisterEntry(__T14, ((void ( *)(const double *, const double *, const int, const bool, const int, const double *, const double *, const double *, double *, double *, const double *, const double *, const double *, const int, const double, const double, const int))gpu_linear_interpolate_spherical_density), _Z40gpu_linear_interpolate_spherical_densityPKdS0_ibiS0_S0_S0_PdS1_S0_S0_S0_iddi, (-1));
__cudaRegisterEntry(__T14, ((void ( *)(const int, const double *, const double *, const double *, const double *, const double *, const double *, const double *, const int *, const int, double *, double *, double *, double *, double *, const int))gpu_make_grid), _Z13gpu_make_gridiPKdS0_S0_S0_S0_S0_S0_PKiiPdS3_S3_S3_S3_i, (-1));
__cudaRegisterVariable(__T14, __shadow_var(gpu_nex,::gpu_nex), 0, 4Ui64, 1, 0);
__cudaRegisterVariable(__T14, __shadow_var(gpu_nmo,::gpu_nmo), 0, 4Ui64, 1, 0);
__cudaRegisterVariable(__T14, __shadow_var(gpu_ncen,::gpu_ncen), 0, 4Ui64, 1, 0);
__cudaRegisterVariable(__T14, __shadow_var(gpu_MaxGrid,::gpu_MaxGrid), 0, 4Ui64, 1, 0);
__cudaRegisterVariable(__T14, __shadow_var(gpu_bragg_angstrom,::gpu_bragg_angstrom), 0, 912Ui64, 1, 0);
__cudaRegisterVariable(__T14, __shadow_var(gpu_type_vector,::gpu_type_vector), 0, 672Ui64, 1, 0);
}
static void __sti____cudaRegisterAll(void)
{
____cudaRegisterLinkedBinary(__nv_cudaEntityRegisterCallback);
}
