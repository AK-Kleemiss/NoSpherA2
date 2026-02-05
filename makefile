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

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

all: NoSpherA2

# --- Cross-platform rust check ---
check_rust:
	@command -v rustc >/dev/null 2>&1 || { \
	  echo "Rust is not installed or not on PATH. Install via https://www.rust-lang.org/tools/install"; exit 1; }; \
	echo "Found: $$(rustc --version)"


ifeq ($(NAME),MAC)
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
featomic_x86_64:
	@$(MAKE) NATIVE_ARCH=x86_64 featomic

featomic_arm64: check_rust
	@$(MAKE) NATIVE_ARCH=arm64 featomic
endif

ifeq ($(NAME),LINUX)
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



define DETECT_MKL
ifeq ($(NAME),LINUX)
  _MKL_MATCHES := $(wildcard \
    $(MKLROOT)/include/mkl.h \
    /opt/intel/oneapi/mkl/latest/include/mkl.h \
    /usr/local/intel/oneapi/mkl/latest/include/mkl.h \
    /usr/intel/oneapi/mkl/latest/include/mkl.h \
  )
  ifneq ($$(_MKL_MATCHES),)
    MKLROOT := $$(patsubst %/include/mkl.h,%,$$(_MKL_MATCHES))
  endif
endif
endef
$(eval $(DETECT_MKL))


ifeq ($(NAME),MAC)
IntelMKL:
	@brew list --versions libomp >/dev/null 2>&1 || \
	  { echo "Installing libomp (arm64)…"; brew install libomp; }
endif

ifeq ($(NAME),LINUX)
IntelMKL:
ifeq ($(strip $(MKLROOT)),)
	@echo MKL not found, building/installing Intel MKL for Linux
	@wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/6a17080f-f0de-41b9-b587-52f92512c59a/intel-onemkl-2025.3.1.11_offline.sh
	@echo Installing MKL, this will take some time! DO NOT CLOSE THE TERMINAL!
	@sh intel-onemkl-2025.3.1.11_offline.sh -a -s --eula accept
	@$(eval $(DETECT_MKL))
else
	@echo Skipping IntelMKL build, found MKL at: $(MKLROOT)
endif
endif



ifeq ($(NAME),MAC)
LibCint:
	@if [ ! -f Lib/LibCint/lib_$(NATIVE_ARCH)/libcint.a ]; then \
		echo 'Building LibCint for $(NATIVE_ARCH), since Lib/LibCint_$(NATIVE_ARCH)/lib/cint.a doesnt exist'; \
		cd libcint && mkdir -p build_$(NATIVE_ARCH) && cd build_$(NATIVE_ARCH) &&\
		cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=../../Lib/LibCint_$(NATIVE_ARCH) -DCMAKE_OSX_ARCHITECTURES=$(NATIVE_ARCH) -DCMAKE_OSX_DEPLOYMENT_TARGET=13.3 .. && \
		make install; \
	else \
		echo 'Skipping LibCint build, Lib/LibCint/lib_$(NATIVE_ARCH)/cint.a already exists'; \
	fi
	
LibCint_x86_64:
	@$(MAKE) NATIVE_ARCH=x86_64 LibCint

LibCint_arm64:
	@$(MAKE) NATIVE_ARCH=arm64 LibCint
endif



ifeq ($(NAME),LINUX)
LibCint:
	@if [ ! -f Lib/LibCint/lib/libcint.a ]; then \
		cd libcint && mkdir -p build && cd build && \
		cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=../../Lib/LibCint .. && \
		make install; \
	else \
		echo 'Skipping LibCint build, Lib\LibCint\lib\cint.a already exists'; \
	fi
endif

ifeq ($(NAME),WINDOWS)
WIN_SLN  ?= Windows/NoSpherA2.sln
WIN_CFG  ?= Release
WIN_PLAT ?= x64
NoSpherA2: 
	echo Building $(WIN_SLN) ($(WIN_CFG) $(WIN_PLAT)) 
	msbuild $(WIN_SLN) /p:Configuration=$(WIN_CFG) /p:Platform=$(WIN_PLAT)

NoSpherA2_Debug:
	@echo Building NoSpherA2_Debug for $(NAME)
	@$(MAKE) WIN_CFG=Debug NoSpherA2

clean:
	@msbuild $(WIN_SLN) /t:Clean /p:Configuration=$(WIN_CFG) /p:Platform=$(WIN_PLAT)
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
