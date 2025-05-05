ifeq ($(OS), Windows_NT)
  NAME := WINDOWS
else
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S), Linux)
    NAME := LINUX
  endif
  ifeq ($(UNAME_S), Darwin)
    NAME := MAC
  endif
endif

all: check_rust NoSpherA2

# Check for Rust
check_rust:
	@rustc --version >nul 2>&1 || (echo Rust is not installed. Please install Rust from https://www.rust-lang.org/tools/install && exit 1)

OpenBLAS:
ifeq ($(NAME),WINDOWS)
	@if not exist OpenBLAS/build/lib/RELEASE/openblas.lib ( \
		@echo Building OpenBLAS for $(NAME) && \
		cd OpenBLAS && mkdir build && cd build && cmake -G "Visual Studio 17 2022" -DCMAKE_BUILD_TYPE=Release -DNOFORTRAN=ON .. && msbuild -nologo OpenBLAS.sln -p:Configuration=Release -m && cmake -DBUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=install -P cmake_install.cmake \
	) else (\
		@echo OpenBLAS already built \
	)
else
	@echo OpenBLAS built handled internally for $(NAME)
endif

featomic: check_rust
ifeq ($(NAME),WINDOWS)
    @if not exist featomic\featomic_install\lib\metatensor.lib ( \
		@echo Building featomic for $(NAME) && \
		cd featomic\featomic && if not exist build mkdir build && cd build && cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../featomic_install --fresh .. && make install \
	) else ( \
		@echo featomic already built \
	)
endif
ifeq ($(NAME),MAC)
	MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
	@if [ ! -f featomic/featomic_install_arm/lib/libfeatomic.a ]; then \
		echo 'Building featomic, since featomic/featomic_install_arm/lib/libfeatomic.a doesnt exist'; \
		cd $(MAKEFILE_DIR)/featomic/featomic && mkdir -p build_arm && cd build_arm && cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON  -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../featomic_install_arm .. && make install\
	else \
		echo 'Skipping featomic build, featomic/featomic_install_arm/lib/libfeatomic.a already exists'; \
	fi
else
	@if [ ! -f featomic/featomic_install/lib/libfeatomic.a ]; then \
		echo 'Building featomic, since featomic/featomic_install/lib/libfeatomic.a doesnt exist'; \
		cd cd featomic/featomic && mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON  -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../featomic_install .. && make install \
	else \
		echo 'Skipping featomic build, featomic/featomic_install/lib/libfeatomic.a already exists'; \
	fi
endif

featomic_x86: check_rust
ifeq ($(NAME),MAC)
	MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
	@if [ ! -f featomic/featomic_install_x86/lib/libfeatomic.a ]; then \
		echo 'Building featomic, since featomic/featomic_install_x86/lib/libfeatomic.a doesnt exist'; \
		rustup target add x86_64-apple-darwin; \
		cd $(MAKEFILE_DIR)/featomic/featomic && mkdir -p build_x86_64 && cd build_x86_64 && cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON  -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../featomic_install_x86 -DCMAKE_OSX_ARCHITECTURES=x86_64 -DRUST_BUILD_TARGET="x86_64-apple-darwin" .. && make install; \
	else \
		echo 'Skipping featomic build, featomic/featomic_install_x86/lib/libfeatomic.a already exists'; \
	fi
else

NoSpherA2_Debug: featomic featomic_x86 OpenBLAS
	@echo Building NoSpherA2_Debug for $(NAME)
ifeq ($(NAME),WINDOWS)
	cd Windows && msbuild NoSpherA2.sln /p:Configuration=Debug /p:Platform=x64 && cd .. && copy Windows\x64\Debug\NoSpherA2.exe .
endif
ifeq ($(NAME),LINUX)
	@echo Start making Linux executable for NoSpherA2_Debug
	rm -f NoSpherA2_Debug
	cd Linux && rm -f NoSpherA2_Debug && make NoSpherA2_Debug -j
endif
ifeq ($(NAME),MAC)
	@echo Start making Mac executable for NoSpherA2_Debug
	cd Mac && rm -f NoSpherA2_Debug && make NoSpherA2_Debug -j
endif

NoSpherA2: featomic featomic_x86 OpenBLAS
	@echo Building NoSpherA2 for $(NAME)
ifeq ($(NAME),WINDOWS)
	cd Windows && msbuild NoSpherA2.sln /p:Configuration=Release /p:Platform=x64 && cd .. && copy Windows\x64\Release\NoSpherA2.exe .
endif
ifeq ($(NAME),LINUX)
	@echo Start making Linux executable
	@rm -f NoSpherA2
	@cd Linux && rm -f NoSpherA2 && make all -j
endif
ifeq ($(NAME),MAC)
	@echo Start making Mac executable
	@rm -f NoSpherA2
	@cd Mac && rm -f NoSpherA2 && make all -j
endif

NoSpherA2_arm: featomic OpenBLAS
ifeq ($(NAME),MAC)
	@echo Start making Mac executable
	@rm -f NoSpherA2
	@cd Mac && rm -f NoSpherA2 && make NoSpherA2_arm -j && cp NoSpherA2_arm ../NoSpherA2
endif

test: NoSpherA2
	make -C tests all -k -B
tests: NoSpherA2
	make -C tests all -k -B


.PHONY: test tests NoSpherA2 all NoSpherA2_Debug OpenBLAS
