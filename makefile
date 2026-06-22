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
export MACOSX_DEPLOYMENT_TARGET=15.0
export CMAKE_OSX_DEPLOYMENT_TARGET=15.0
export RUSTFLAGS=-C link-arg=-mmacosx-version-min=15.0
endif

BUILD_TESTS ?= 0
CONFIG ?= Release

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
		make install && cd $(MAKEFILE_DIR)/ && \
		if [ -f Lib/featomic_install/lib64/libmetatensor.a ]; then \
			echo "Copying lib64/libmetatensor.a to lib/"; \
		    cp Lib/featomic_install/lib64/libmetatensor.a Lib/featomic_install/lib; \
		fi; \
	else \
		echo 'Skipping featomic build, Lib/featomic_install/lib/libfeatomic.a already exists'; \
	fi

endif



ifeq ($(NAME),MAC)
IntelMKL:
	@brew list --versions libomp >/dev/null 2>&1 || \
	  { echo "Installing libomp (arm64)…"; brew install libomp; }
endif

ifeq ($(NAME),LINUX)
MKLROOT := $(CURDIR)/Lib/MKL
IntelMKL:
ifeq ($(wildcard $(MKLROOT)),)
	@echo MKL not found, building/installing Intel MKL for Linux
	@wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/6a17080f-f0de-41b9-b587-52f92512c59a/intel-onemkl-2025.3.1.11_offline.sh
	@echo Installing MKL, this will take some time! DO NOT CLOSE THE TERMINAL!
	@sh intel-onemkl-2025.3.1.11_offline.sh -a -s --install-dir=$(MKLROOT) --eula accept
else
	@echo Skipping IntelMKL build, found MKL at: $(MKLROOT)
endif
endif

BasisSetConverter: 
	@if [ ! -f Lib/BasisSetConverter ]; then \
		echo 'Building BasisSetConverter, since Lib/BasisSetConverter doesnt exist'; \
		cd BasisSetGenerator && mkdir -p build \
		&& cmake -S . -B build -DCMAKE_OSX_ARCHITECTURES=$(NATIVE_ARCH)\
    	&& cmake --build build --config Release \
    	&& cmake --install build --config Release --prefix "../Lib"; \
	else \
		echo 'Skipping BasisSetConverter build, found Lib/BasisSetConverter'; \
	fi

ifeq ($(NAME),MAC)
LibCint:
	@if [ ! -f Lib/LibCint_$(NATIVE_ARCH)/lib/libcint.a ]; then \
		echo 'Building LibCint for $(NATIVE_ARCH), since Lib/LibCint_$(NATIVE_ARCH)/lib/libcint.a doesnt exist'; \
		cd libcint && mkdir -p build_$(NATIVE_ARCH) && cd build_$(NATIVE_ARCH) &&\
		cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=../../Lib/LibCint_$(NATIVE_ARCH) \
			-DENABLE_STATIC="ON" -DWITH_CINT2_INTERFACE="OFF" -DBUILD_SHARED_LIBS="OFF" -DPYPZPX="ON" \
			-DWITH_RANGE_COULOMB="ON" -DCMAKE_OSX_ARCHITECTURES=$(NATIVE_ARCH) -DCMAKE_OSX_DEPLOYMENT_TARGET=13.3 .. && \
		make install; \
	else \
		echo 'Skipping LibCint build, Lib/LibCint_$(NATIVE_ARCH)/lib/libcint.a already exists'; \
	fi

occ: LibCint
	@if [ ! -f Lib/occ_$(NATIVE_ARCH)/lib/libocc.a ]; then \
		echo 'Building OCC for $(NATIVE_ARCH), since Lib/occ_$(NATIVE_ARCH)/lib/libocc.a doesnt exist'; \
        cmake --workflow --preset macos-release-$(NATIVE_ARCH) && \
        cmake --install ./build-macos-release-$(NATIVE_ARCH) && \
		libtool -static -o Lib/occ_$(NATIVE_ARCH)/lib/libocc.a Lib/occ_$(NATIVE_ARCH)/lib/*.a; \
	else \
		echo 'Skipping occ build, Lib/occ_$(NATIVE_ARCH)/lib/libocc.a already exists'; \
	fi

LibCint_x86_64:
	@$(MAKE) NATIVE_ARCH=x86_64 LibCint

LibCint_arm64:
	@$(MAKE) NATIVE_ARCH=arm64 LibCint

occ_x86_64: LibCint_x86_64
	@$(MAKE) NATIVE_ARCH=x86_64 occ

occ_arm64: LibCint_arm64
	@$(MAKE) NATIVE_ARCH=arm64 occ

BasisSetConverter_x86_64:
	@$(MAKE) NATIVE_ARCH=x86_64 BasisSetConverter

BasisSetConverter_arm64:
	@$(MAKE) NATIVE_ARCH=arm64 BasisSetConverter
endif



ifeq ($(NAME),LINUX)
LibCint:
	@if [ ! -f Lib/LibCint/lib/libcint.a ]; then \
		cd libcint && mkdir -p build && cd build && \
		cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE \
			-DENABLE_STATIC="ON" -DWITH_CINT2_INTERFACE="OFF" -DBUILD_SHARED_LIBS="OFF" -DPYPZPX="ON" \
            -DWITH_RANGE_COULOMB="ON" -DCMAKE_INSTALL_PREFIX=../../Lib/LibCint .. && \
		make install && cd $(MAKEFILE_DIR)/ && \
		if [ -f Lib/LibCint/lib64/libcint.a ]; then \
			echo "Copying lib64/libcint.a to lib/"; \
			mkdir -p Lib/LibCint/lib && \
		    cp Lib/LibCint/lib64/libcint.a Lib/LibCint/lib; \
		fi; \
	else \
		echo 'Skipping LibCint build, Lib/LibCint/lib/libcint.a already exists'; \
	fi

occ: LibCint
	@if [ ! -f Lib/occ/lib/libocc.a ]; then \
		echo 'Building OCC, since Lib/occ/lib/libocc_main.a doesnt exist'; \
		cmake --workflow --preset linux-occ-gcc && \
		cmake --install build-linux-occ-gcc && \
	  	cd build-linux-occ-gcc && ar M < libocc.ar && cd .. && \
		cp build-linux-occ-gcc/liblibocc.a Lib/occ/lib/libocc.a; \
	else \
		echo 'Skipping OCC build, Lib/occ/lib/libocc_main.a already exists'; \
	fi

endif


ifeq ($(NAME),WINDOWS)
WIN_SLN  ?= Windows/NoSpherA2.sln
WIN_CFG  ?= Release
WIN_PLAT ?= x64
WIN_VCTOOLSVERSION ?= 14.51.36231
NoSpherA2: 
	echo Building $(WIN_SLN) ($(WIN_CFG) $(WIN_PLAT)) 
	msbuild $(WIN_SLN) /p:Configuration=$(WIN_CFG) /p:Platform=$(WIN_PLAT) /p:VCToolsVersion=$(WIN_VCTOOLSVERSION)

NoSpherA2_Debug:
	@echo Building NoSpherA2_Debug for $(NAME)
	@$(MAKE) WIN_CFG=Debug NoSpherA2

clean:
	@msbuild $(WIN_SLN) /t:Clean /p:Configuration=$(WIN_CFG) /p:Platform=$(WIN_PLAT) /p:VCToolsVersion=$(WIN_VCTOOLSVERSION)
endif

ifeq ($(NAME),LINUX)
NoSpherA2: IntelMKL featomic LibCint occ BasisSetConverter
	@echo Building NoSpherA2 for $(NAME) with $(CONFIG) configuration
	@cd Linux && $(MAKE) all -j BUILD_TESTS=$(BUILD_TESTS) CONFIG=$(CONFIG) 

NoSpherA2_Debug: 
	@$(MAKE) BUILD_TESTS=$(BUILD_TESTS) CONFIG=Debug NoSpherA2 

clean:
	@cd Linux && make clean
endif

ifeq ($(NAME),MAC)
MAC_OUT_DIR := build/$(CONFIG)/bin

NoSpherA2: IntelMKL featomic LibCint occ BasisSetConverter
	@echo Building NoSpherA2 for $(NAME) with $(CONFIG) configuration
	@cd Mac && $(MAKE) NoSpherA2 -j CONFIG=$(CONFIG) ARCH=$(NATIVE_ARCH)
	@echo 'Done building NoSpherA2 (see $(MAC_OUT_DIR))'

NoSpherA2_arm64: IntelMKL featomic_arm64 LibCint_arm64 occ_arm64 BasisSetConverter_arm64
	@$(MAKE) NATIVE_ARCH=arm64 NoSpherA2

NoSpherA2_x86_64: IntelMKL featomic_x86_64 LibCint_x86_64 occ_x86_64 BasisSetConverter_x86_64
	@$(MAKE) NATIVE_ARCH=x86_64 NoSpherA2

NoSpherA2_lipo: NoSpherA2_arm64 NoSpherA2_x86_64
	@echo Building universal NoSpherA2 for $(NAME)
	@cd Mac && $(MAKE) lipo -j CONFIG=$(CONFIG)
	@echo 'Done building NoSpherA2 (see $(MAC_OUT_DIR))'

clean:
	@cd Mac && make clean
endif

test: NoSpherA2
	make -C tests all -k -B
tests: NoSpherA2
	make -C tests all -k -B


.PHONY: test tests NoSpherA2 all NoSpherA2_Debug clean IntelMKL featomic check_rust LibCint
