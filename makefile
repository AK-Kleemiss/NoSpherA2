ifeq ($(OS), Windows_NT)
  NAME := WINDOWS
else
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S), Linux)
	COMP := g++
	C_COMP := gcc
	NAME := LINUX
  endif
  ifeq ($(UNAME_S), Darwin)
	COMP := clang++
	C_COMP := clang
	NAME := MAC
	# Detect native architecture on macOS
	NATIVE_ARCH := $(shell uname -m)
	ifeq ($(NATIVE_ARCH),arm64)
		NATIVE_ARCH := arm64
	else
		NATIVE_ARCH := x86_64
	endif
  endif
endif

#Set some environemt variables for macOS builds
ifeq ($(NAME),MAC)
export MACOSX_DEPLOYMENT_TARGET=13.3
export CMAKE_OSX_DEPLOYMENT_TARGET=13.3
export RUSTFLAGS=-C link-arg=-mmacosx-version-min=13.3
endif

ifeq ($(OS),Windows_NT)

# Find latest VS installation path via vswhere (PowerShell used only to query)
VS_INSTALL := $(strip $(shell powershell -NoProfile -Command "& '%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe' -latest -products * -property installationPath" ))

VCVARSALL   := $(VS_INSTALL)\VC\Auxiliary\Build\vcvarsall.bat
#Print found VS path
ifeq ($(VS_INSTALL),)
$(error "Visual Studio not found. Please install Visual Studio 2022 with C++ workload.")
else
$(info Found Visual Studio installation at: $(VS_INSTALL))
endif
endif

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

all: check_rust NoSpherA2

# --- Cross-platform rust check ---
ifeq ($(OS),Windows_NT)
check_rust:
	@powershell -NoProfile -Command "if (Get-Command rustc -ErrorAction SilentlyContinue) { rustc --version } else { Write-Error 'Rust not found install via https://www.rust-lang.org/tools/install'; exit 1 }"
else
check_rust:
	@command -v rustc >/dev/null 2>&1 || { \
	  echo "Rust is not installed or not on PATH. Install via https://www.rust-lang.org/tools/install"; exit 1; }; \
	echo "Found: $$(rustc --version)"
endif

ifeq ($(NAME),WINDOWS)
featomic: check_rust
	@if not exist Lib\featomic_install\lib\metatensor.lib ( \
		echo Building featomic for $(NAME) && \
		cd $(MAKEFILE_DIR)/featomic/featomic && \
		(if exist build rd /s /q build) && \
		mkdir build && \
		cd build && \
		echo Starting build && \
		cmake -G "Visual Studio 17 2022" -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="../../../Lib/featomic_install" .. && \
		msbuild -nologo .\featomic.sln -p:Configuration=Release && \
		msbuild -nologo .\INSTALL.vcxproj -p:Configuration=Release \
	) else ( \
		echo featomic already built \
	)
else ifeq ($(NAME),MAC)
featomic:  check_rust
	@if [ ! -f Lib/featomic_install_$(NATIVE_ARCH)/lib/libfeatomic.a ]; then \
		echo 'Building featomic for $(NATIVE_ARCH), since Lib/featomic_install_$(NATIVE_ARCH)/lib/libfeatomic.a doesnt exist'; \
		cd $(MAKEFILE_DIR)/featomic/featomic && \
		mkdir -p build_$(NATIVE_ARCH) && \
		cd build_$(NATIVE_ARCH) && \
		if [ "$(NATIVE_ARCH)" = "x86_64" ]; then \
			rustup target add x86_64-apple-darwin; \
			cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON -DCMAKE_OSX_ARCHITECTURES=$(NATIVE_ARCH) -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../../Lib/featomic_install_$(NATIVE_ARCH) -DRUST_BUILD_TARGET="x86_64-apple-darwin" .. ; \
		else \
			cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON -DCMAKE_OSX_ARCHITECTURES=$(NATIVE_ARCH) -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../../Lib/featomic_install_$(NATIVE_ARCH) .. ; \
		fi && \
		make install || true; \
	else \
		echo 'Skipping featomic build, Lib/featomic_install_$(NATIVE_ARCH)/lib/libfeatomic.a already exists'; \
	fi
else
featomic: check_rust
	@if [ ! -f Lib/featomic_install/lib/libfeatomic.a ]; then \
		echo 'Building featomic, since Lib/featomic_install/lib/libfeatomic.a doesnt exist'; \
		cd $(MAKEFILE_DIR)/featomic/featomic && \
		mkdir -p build && \
		cd build && \
		cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON  -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../../Lib/featomic_install .. && \
		make install; \
	else \
		echo 'Skipping featomic build, Lib/featomic_install/lib/libfeatomic.a already exists'; \
	fi
endif

featomic_x86_64:
	@$(MAKE) NATIVE_ARCH=x86_64 featomic

featomic_arm64: check_rust
	@$(MAKE) NATIVE_ARCH=arm64 featomic


intel_ROOT := $(CURDIR)/Lib/MKL
IntelMKL:
ifeq ($(NAME),WINDOWS)
	MSBuild .\Windows\NoSpherA2.sln /t:Restore /p:RestorePackagesPath=Windows/packages /p:RestoreTimeout=300000
else ifeq ($(NAME),MAC)
	@brew list --versions libomp >/dev/null 2>&1 || \
	  { echo "Installing libomp (arm64)…"; brew install libomp; }
else
	@if [ ! -f Lib/MKL/mkl/2025.2/lib/libmkl_core.a ]; then \
		wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/47c7d946-fca1-441a-b0df-b094e3f045ea/intel-onemkl-2025.2.0.629_offline.sh; \
		sh intel-onemkl-2025.2.0.629_offline.sh -a -s --install-dir=$(intel_ROOT) --eula accept; \
	else \
		echo 'Skipping IntelMKL build, MKL\\mkl\\2025.2\\lib\\libmkl_core.a already exists'; \
	fi
endif


ifeq ($(NAME),WINDOWS)
LibCint:
	@if not exist Lib\LibCint\lib\cint.lib ( \
		echo Building LibCint for $(NAME) && \
		@cd libcint && (if not exist build mkdir build) && cd build && \
		cmd /S /C "call "$(VCVARSALL)" x64 >nul && cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=../../Lib/LibCint -DCMAKE_C_COMPILER=cl .." && \
		cmake --build . --config RELEASE && cmake --install . \
	) else ( \
		echo Skipping LibCint build, Lib\LibCint\lib\cint.lib already exists \
	)
else ifeq ($(NAME),MAC)
LibCint:
	@if [ ! -f Lib/LibCint/lib_$(NATIVE_ARCH)/cint.a ]; then \
		echo 'Building LibCint for $(NATIVE_ARCH), since Lib/LibCint_$(NATIVE_ARCH)/lib/cint.a doesnt exist'; \
		cd libcint && mkdir -p build_$(NATIVE_ARCH) && cd build_$(NATIVE_ARCH) &&\
		cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=../../Lib/LibCint -DCMAKE_OSX_ARCHITECTURES=$(NATIVE_ARCH) -DCMAKE_OSX_DEPLOYMENT_TARGET=13.3 .. && \
		make install; \
	else \
		echo 'Skipping LibCint build, Lib/LibCint/lib_$(NATIVE_ARCH)/cint.a already exists'; \
	fi
else
LibCint:
	@if [ ! -f Lib/LibCint/lib/cint.a ]; then \
		cd libcint && mkdir -p build && cd build && \
		cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=../../Lib/LibCint .. && \
		make install; \
	else \
		echo 'Skipping LibCint build, Lib\LibCint\lib\cint.a already exists'; \
	fi
endif

LibCint_x86_64:
	@$(MAKE) NATIVE_ARCH=x86_64 LibCint

LibCint_arm64:
	@$(MAKE) NATIVE_ARCH=arm64 LibCint

ifeq ($(NAME),WINDOWS)
NoSpherA2: IntelMKL featomic LibCint
	@cd Windows && msbuild NoSpherA2.sln /p:Configuration=Release /p:Platform=x64 && cd .. && copy Windows\x64\Release\NoSpherA2.exe . && copy Windows\x64\Release\libiomp5md.dll

NoSpherA2_Debug: IntelMKL featomic LibCint
	@echo Building NoSpherA2_Debug for $(NAME)
	@cd Windows && msbuild NoSpherA2.sln /p:Configuration=Debug /p:Platform=x64 && cd .. && copy Windows\x64\Debug\NoSpherA2.exe .

clean:
	@cd Windows && if exist x64 rd /s /q x64
endif

ifeq ($(NAME),LINUX)
NoSpherA2: IntelMKL featomic LibCint
	@echo Start making Linux executable
	@rm -f NoSpherA2
	@cd Linux && rm -f NoSpherA2 && make all -j

NoSpherA2_Debug: IntelMKL featomic LibCint
	@echo Building NoSpherA2_Debug for $(NAME)
	@rm -f NoSpherA2_Debug
	@cd Linux && rm -f NoSpherA2_Debug && make NoSpherA2_Debug -j

clean:
	@cd Linux && make clean

endif

ifeq ($(NAME),MAC)
NoSpherA2: IntelMKL featomic LibCint
	@echo Start making Mac $(NATIVE_ARCH) executable
	@rm -f NoSpherA2_$(NATIVE_ARCH)
	@cd Mac && rm -f NoSpherA2_$(NATIVE_ARCH) && make NoSpherA2_$(NATIVE_ARCH) -j && cp NoSpherA2_$(NATIVE_ARCH) ../NoSpherA2

NoSpherA2_arm64: IntelMKL featomic_arm64 LibCint_arm64
	@echo Start making Mac arm64 executable
	@rm -f NoSpherA2_arm64
	@cd Mac && rm -f NoSpherA2_arm64 && make NoSpherA2_arm64 -j && cp NoSpherA2_arm64 ../NoSpherA2

NoSpherA2_x86_64: IntelMKL featomic_x86_64 LibCint_x86_64
	@echo Start making Mac x86_64 executable
	@rm -f NoSpherA2_x86_64
	@cd Mac && rm -f NoSpherA2_x86_64 && make NoSpherA2_x86_64 -j && cp NoSpherA2_x86_64 ../NoSpherA2_x86_64

NoSpherA2_lipo: IntelMKL featomic_arm64 featomic_x86_64 LibCint_arm64 LibCint_x86_64
	@echo Start making Mac universal executable
	@rm -f NoSpherA2
	@cd Mac && rm -f NoSpherA2 && make NoSpherA2 -j

clean:
	@cd Mac && make clean
endif

test: NoSpherA2
	make -C tests all -k -B
tests: NoSpherA2
	make -C tests all -k -B


.PHONY: test tests NoSpherA2 all NoSpherA2_Debug clean IntelMKL featomic check_rust LibCint
