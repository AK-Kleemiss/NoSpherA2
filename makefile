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
  endif
endif

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

all: check_rust NoSpherA2

# Check for Rust
check_rust:
	@rustc --version >nul 2>&1 || (echo Rust is not installed. Please install Rust from https://www.rust-lang.org/tools/install && exit 1)

featomic: check_rust
ifeq ($(NAME),WINDOWS)
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
	@if [ ! -f Lib/featomic_install/lib/libfeatomic.a ]; then \
		echo 'Building featomic, since Lib/featomic_install/lib/libfeatomic.a doesnt exist'; \
		cd $(MAKEFILE_DIR)/featomic/featomic && \
		mkdir -p build && \
		cd build && \
		cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON  -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../../Lib/featomic_install .. && \
		make install || true;
	else \
		echo 'Skipping featomic build, Lib/featomic_install/lib/libfeatomic.a already exists'; \
	fi
else
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

featomic_x86: check_rust
ifeq ($(NAME),MAC)
	MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
	@if [ ! -f Lib/featomic_install_x86/lib/libfeatomic.a ]; then \
		echo 'Building featomic, since Lib/featomic_install_x86/lib/libfeatomic.a doesnt exist'; \
		rustup target add x86_64-apple-darwin; \
		cd $(MAKEFILE_DIR)/featomic/featomic && \
		mkdir -p build_x86_64 && \
		cd build_x86_64 && \
		cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=../../../Lib/featomic_install_x86 -DCMAKE_OSX_ARCHITECTURES=x86_64 -DRUST_BUILD_TARGET="x86_64-apple-darwin" .. && \
		make install; \
	else \
		echo 'Skipping featomic build, featomic/featomic_install_x86/lib/libfeatomic.a already exists'; \
	fi
else
	@echo No need for cross compile featomic on non-MacOS systems
endif

IntelMKL:
ifeq ($(NAME),WINDOWS)
	MSBuild .\Windows\NoSpherA2.sln /t:Restore -p:RestorePackagesConfig=true
else ifeq ($(NAME),MAC)
	echo 'Not implemented'
else
	echo 'Not implemented'
endif


NoSpherA2_Debug: IntelMKL featomic featomic_x86
	@echo Building NoSpherA2_Debug for $(NAME)
ifeq ($(NAME),WINDOWS)
	@cd Windows && msbuild NoSpherA2.sln /p:Configuration=Debug /p:Platform=x64 && cd .. && copy Windows\x64\Debug\NoSpherA2.exe .
endif
ifeq ($(NAME),LINUX)
	@echo Start making Linux executable for NoSpherA2_Debug
	@rm -f NoSpherA2_Debug
	@cd Linux && rm -f NoSpherA2_Debug && make NoSpherA2_Debug -j
endif
ifeq ($(NAME),MAC)
	@echo Start making Mac executable for NoSpherA2_Debug
	@cd Mac && rm -f NoSpherA2_Debug && make NoSpherA2_Debug -j
endif

NoSpherA2: IntelMKL featomic
	@echo Building NoSpherA2 for $(NAME)
ifeq ($(NAME),WINDOWS)
	@cd Windows && msbuild NoSpherA2.sln /p:Configuration=Release /p:Platform=x64 && cd .. && copy Windows\x64\Release\NoSpherA2.exe .
endif
ifeq ($(NAME),LINUX)
	@echo Start making Linux executable
	@rm -f NoSpherA2
	@cd Linux && rm -f NoSpherA2 && make all -j
endif
ifeq ($(NAME),MAC)
	@echo Start making Mac executable
	@rm -f NoSpherA2
	@cd Mac && rm -f NoSpherA2_native && make NoSpherA2_native -j
endif

NoSpherA2_arm: IntelMKL featomic
ifeq ($(NAME),MAC)
	@echo Start making Mac executable
	@rm -f NoSpherA2
	@cd Mac && rm -f NoSpherA2 && make NoSpherA2_arm -j && cp NoSpherA2_arm ../NoSpherA2
endif

NoSpherA2_x86: IntelMKL featomic_x86
ifeq ($(NAME),MAC)
	@echo Start making Mac executable
	@rm -f NoSpherA2
	@cd Mac && rm -f NoSpherA2 && make NoSpherA2_x86 -j && cp NoSpherA2_x86 ../NoSpherA2
endif

NoSpherA2_lipo:IntelMKL featomic featomic_x86
ifeq ($(NAME),MAC)
	@echo Start making Mac executable
	@rm -f NoSpherA2
	@cd Mac && rm -f NoSpherA2 && make NoSpherA2 -j && cp NoSpherA2 ../NoSpherA2
endif

test: NoSpherA2
	make -C tests all -k -B
tests: NoSpherA2
	make -C tests all -k -B


.PHONY: test tests NoSpherA2 all NoSpherA2_Debug
