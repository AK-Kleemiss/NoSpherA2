# Updates to NoSpherA2:
## Added OCC Support
### 1. Added new way to download MKL using conda. 
It just downloads the necessary packages, avoiding the full intel installation.
It uses micromamba (a much faster version of miniconda written in C++) to create an environment and install the packages.

This is also good because it works to install multiple packages and it is possible to change versions of MKL easily.
Also, it is possible to install other packages. It makes caching easy in Actions because you only have to hash the 
environment.yml file.

### 2. Unified dealing with dependencies
All dependencies that are to be built from source are available in the 3rdparty directory. This also simplifies actions
as the 3rdparty/CMakeLists.txt can be hashed for caching.

### 3. Refactored actions to use CMake + Ninja
The actions were refactored to use CMake + Ninja. This makes it easier to add new configurations and it is faster than
MSbuild. The mamba environment is also cached. Also, sccache is saved after each build. Some things still need to be improved, 
like caching the 3rdparty builds, but it is a good start. The benefit of this is that you can automatically invalidate the
caches by hashing only the 3rdpary CMakeLists.txt file. That way, the dependencies can be rebuilt on updates and tested
more robustly.

### 3. Added support for OCC
Support for OCC was added. It is possible to build OCC using 
```
make.exe occ
```
in the root directory. I wasn't able to make it work well with NoSpherA2 vcxproj and sln files, so if someone could help it
would be appreciated.

It is still possible to build NoSpherA2 using CMake and Ninja:
```
choco install ninja sccache
```
```
cmake .. -G "Ninja"
    -DCMAKE_BUILD_TYPE=Release `
    -DCMAKE_C_COMPILER=clang-cl.exe `
    -DCMAKE_LINKER=lld-link `
    -DCMAKE_CXX_COMPILER=clang-cl.exe `
    -DCMAKE_EXE_LINKER_FLAGS="runtimeobject.lib" `
    -DCMAKE_SHARED_LINKER_FLAGS="runtimeobject.lib" `
    -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded `
    -DBUILD_SHARED_LIBS=OFF `
    -DTBB_BUILD_SHARED=OFF
```
This also supports the flags:
```
-DCMAKE_C_COMPILER_LAUNCHER=sccache
-DCMAKE_CXX_COMPILER_LAUNCHER=sccache
```
Sccache is a compiler cache that speeds up recompilation. It runs after preprocessing and stores a hashmap between files
and compiled objects. It is very useful when switching branches often. It also knows when you are building debug and release.

In other platforms the build should work with cmake. Just copying and pasting the command from actions should work.

### 4. Added OCC wavefunction calculation
OCC wavefunction calculations works by creating an OCC input and using the flag -occ input.toml. It forwards the flag
to the normal SCF routine and runs as normal. It translates directly between objects in less than 100ms. It is still kind
of slow and it would be possible to paralelize the for loops, but I thought it was good enough for now. 

The molecular orbital coefficients are written in OCC as a linear array. Probably to avoid the padding in the memory 
that comes from arrays of structs/classes.

A block has to be selected based on the offset of the current molecular orbital. The contraction coefficients from occ are
normalized, as is done in OCC before it outputs to FCHK and then, as it is done in the FCHK reading part of the NoSpherA2 code,
for each shell, the following is done:
$$
MOc_\text{NOS}=\text{vec}\left(|C\rangle \langle \text{MO}c_\text{OCC}|A^\dagger\right)
$$
Where $|C\rangle$ are the contraction coefficients from OCC, $\langle \text{MO}c_\text{OCC}|$ is the block of molecular
orbital coefficients from OCC and $A^\dagger$ is the transformation matrix from spherical to cartesian harmonics transposed.
vec is the operation of turning a matrix into a vector by stacking its columns.

As Eigen stores everything in Column-Major, it is as simple as directly writing to the memory. The size assertions in
Eigen will guarantee that the size is correct. Eigen also makes it easier by offering a Map to do that:
```c++
for (int n=0; n<occ_WF.nbf; ++n)
{
    // ...
    MOs[n].assign_coefficients_size(nex);
    coeffs_ptr = MOs[n].get_coefficient_ptr();
    for (const auto & shell : shells)
    {
        // ...
        const VectorXd& contraction_coeffs = scalar*(2*exp_arr).pow(p)*CC;
        const auto& MOc = mo_go.C.block(basis_offset, n, n_sph, 1);
        Map<MatrixXd> dest_block(coeffs_ptr + write_cursor, n_prim, n_cart);
        dest_block = contraction_coeffs*(A * MOc).transpose();
       // ... 
    }
// ...
```
Most other things can be precalculated. This was mainly done for clarity and performance, as it reduces a lot the amount
of nested for loops. As the matrixes aren't that big, each line of them is reduced to at most 3 vector instructions for
multiplication (Using AVX with double precision). AVX2 would reduce most of them to 1 vector instruction. Assuming 
Eigen is working well. I don't think parallelizing this kind of matrix operations is worth it, considering the overhead 
of threads. But the outer loops could be parallelized, I just didn't think it was worth it for now.

Some things can be easily improved, like caching the scalar part of the normalization. But as it was already on par with
reading FCHK using NoSpherA2, I didn't think it was worth it for now.