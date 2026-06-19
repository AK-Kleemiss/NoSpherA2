include common.mk

# High-level build entrypoint for Linux builds.
#
# Targets:
#   - static       : build core static library (build/<CONFIG>/lib/libnosphera2.a)
#   - shared       : build shared object (build/<CONFIG>/bin/libNoSpherA2.so)
#   - NoSpherA2    : build CLI executable (build/<CONFIG>/bin/NoSpherA2)
#   - all          : build static, shared, then executable

REPO_ROOT ?= ..
CONFIG ?= Release
BUILD_DIR ?= $(REPO_ROOT)/build/$(CONFIG)

STATIC_DIR := NoSpherA2_STATIC
EXE_DIR := NoSpherA2
SHARED_DIR := NoSpherA2_SHARED
TEST_DIR := NoSpherA2_Tests
LIB_NOS := ../build/$(CONFIG)/lib/libnosphera2.a

all: NoSpherA2

ifeq ($(BUILD_TESTS),1)
all: tests
endif

static: $(LIB_NOS)

$(LIB_NOS):
	@$(MAKE) -C $(STATIC_DIR) REPO_ROOT=../.. CONFIG=$(CONFIG) all -j


NoSpherA2: static
	@$(MAKE) -C $(EXE_DIR) REPO_ROOT=../.. CONFIG=$(CONFIG) -j
	@mkdir -p ../build/$(CONFIG)/bin
	@cp -f ../Lib/MKL/compiler/latest/lib/libiomp5.so ../build/$(CONFIG)/bin/ 2>/dev/null || true
	@cp -fP ../Lib/occ/lib/libtbb.so* ../build/$(CONFIG)/bin/ 2>/dev/null || true
	@echo 'Done building NoSpherA2 (see build/$(CONFIG)/bin)'

tests: static
	@$(MAKE) -C $(TEST_DIR) REPO_ROOT=../.. CONFIG=$(CONFIG) all -j

clean:
	@$(MAKE) -C $(STATIC_DIR) clean
	@$(MAKE) -C $(EXE_DIR) clean
	@$(MAKE) -C $(TEST_DIR) clean
	-@echo ' '

.PHONY: all static tests clean NoSpherA2 NoSpherA2_Debug



