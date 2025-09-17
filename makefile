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

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

all: check_rust NoSpherA2

# --- Cross-platform rust check ---
ifeq ($(OS),Windows_NT)
check_rust:
	@where rustc >NUL 2>&1 && ( \
	  for /f "usebackq delims=" %%V in (`rustc --version`) do @echo Found: %%V \
	) || ( \
	  echo Rust is not installed or not on PATH. Install via https://www.rust-lang.org/tools/install & exit /b 1 \
	)
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
featomic: check_rust featomic_$(NATIVE_ARCH)
	@echo Built featomic for $(NATIVE_ARCH)
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

featomic_x86_64: check_rust
	@if [ ! -f Lib/featomic_install_x86/lib/libfeatomic.a ]; then \
		echo 'Building featomic for x86_64, since Lib/featomic_install_x86/lib/libfeatomic.a doesnt exist'; \
		rustup target add x86_64-apple-darwin; \
		cd $(MAKEFILE_DIR)/featomic/featomic && \
		mkdir -p build_x86_64 && \
		cd build_x86_64 && \
		cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../../Lib/featomic_install_x86 -DCMAKE_OSX_ARCHITECTURES=x86_64 -DRUST_BUILD_TARGET="x86_64-apple-darwin" .. && \
		make install; \
	else \
		echo 'Skipping featomic build, Lib/featomic_install_x86/lib/libfeatomic.a already exists'; \
	fi

featomic_arm64: check_rust
	@if [ ! -f Lib/featomic_install/lib/libfeatomic.a ]; then \
		echo 'Building featomic, since Lib/featomic_install/lib/libfeatomic.a doesnt exist'; \
		cd $(MAKEFILE_DIR)/featomic/featomic && \
		mkdir -p build && \
		cd build && \
		cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON  -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../../Lib/featomic_install .. && \
		make install || true; \
	else \
		echo 'Skipping featomic build, Lib/featomic_install/lib/libfeatomic.a already exists'; \
	fi

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
		echo 'Skipping IntelMKL build, MKL\mkl\2025.2\lib\libmkl_core.a already exists'; \
	fi
endif


ifeq ($(NAME),WINDOWS)
NoSpherA2: IntelMKL featomic
	@cd Windows && msbuild NoSpherA2.sln /p:Configuration=Release /p:Platform=x64 && cd .. && copy Windows\x64\Release\NoSpherA2.exe . && copy Windows\x64\Release\libiomp5md.dll

NoSpherA2_Debug: IntelMKL featomic
	@echo Building NoSpherA2_Debug for $(NAME)
	@cd Windows && msbuild NoSpherA2.sln /p:Configuration=Debug /p:Platform=x64 && cd .. && copy Windows\x64\Debug\NoSpherA2.exe .

clean:
	@cd Windows && if exist x64 rd /s /q x64
endif

ifeq ($(NAME),LINUX)
NoSpherA2: IntelMKL featomic
	@echo Start making Linux executable
	@rm -f NoSpherA2
	@cd Linux && rm -f NoSpherA2 && make all -j

NoSpherA2_Debug: IntelMKL featomic
	@echo Building NoSpherA2_Debug for $(NAME)
	@rm -f NoSpherA2_Debug
	@cd Linux && rm -f NoSpherA2_Debug && make NoSpherA2_Debug -j

clean:
	@cd Linux && make clean

endif

ifeq ($(NAME),MAC)
NoSpherA2: IntelMKL featomic_$(NATIVE_ARCH)
	@echo Start making Mac $(NATIVE_ARCH) executable
	@rm -f NoSpherA2_$(NATIVE_ARCH)
	@cd Mac && rm -f NoSpherA2_$(NATIVE_ARCH) && make NoSpherA2_$(NATIVE_ARCH) -j && cp NoSpherA2_$(NATIVE_ARCH) ../NoSpherA2

NoSpherA2_arm64: IntelMKL featomic_arm64
	@echo Start making Mac arm64 executable
	@rm -f NoSpherA2_arm64
	@cd Mac && rm -f NoSpherA2_arm64 && make NoSpherA2_arm64 -j && cp NoSpherA2_arm64 ../NoSpherA2

NoSpherA2_x86_64: IntelMKL featomic_x86_64
	@echo Start making Mac x86_64 executable
	@rm -f NoSpherA2_x86_64
	@cd Mac && rm -f NoSpherA2_x86_64 && make NoSpherA2_x86_64 -j && cp NoSpherA2_x86_64 ../NoSpherA2_x86_64

NoSpherA2_lipo: IntelMKL featomic_arm64 featomic_x86_64
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


.PHONY: test tests NoSpherA2 all NoSpherA2_Debug clean IntelMKL featomic check_rust
