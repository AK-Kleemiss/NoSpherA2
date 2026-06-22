include common.mk

# High-level build entrypoint for macOS builds.
# Mirrors Linux/makefile structure:
#   - static : build core static library
#   - NoSpherA2 : build CLI executable
#   - tests : build native test binary (if present/usable)
#   - lipo : build universal binary (arm64 + x86_64)

REPO_ROOT ?= ..
CONFIG ?= Release

BUILD_DIR ?= $(REPO_ROOT)/build/$(CONFIG)

STATIC_DIR := NoSpherA2_STATIC
EXE_DIR := NoSpherA2
TEST_DIR := NoSpherA2_Tests

all: NoSpherA2

# ---- Per-arch builds ----
static:
	@$(MAKE) -C $(STATIC_DIR) REPO_ROOT=../.. CONFIG=$(CONFIG) NATIVE_ARCH=$(ARCH) all -j

NoSpherA2: static
	@$(MAKE) -C $(EXE_DIR) REPO_ROOT=../.. CONFIG=$(CONFIG) NATIVE_ARCH=$(ARCH) all -j
	@echo 'Done building NoSpherA2 (see build/$(CONFIG)/bin)'

tests: static
	@$(MAKE) -C $(TEST_DIR) REPO_ROOT=../.. CONFIG=$(CONFIG) NATIVE_ARCH=$(ARCH) all -j

# ---- Universal binary (arm64 + x86_64) ----
# This produces build/<CONFIG>/bin/NoSpherA2 as a universal binary.
# It expects both per-arch outputs to exist.

lipo: NoSpherA2_arm64 NoSpherA2_x86_64
	@echo 'Lipo → build/$(CONFIG)/bin/NoSpherA2'
	@mkdir -p "$(BUILD_DIR)/bin"
	@lipo -create -output "$(BUILD_DIR)/bin/NoSpherA2" \
	  "$(BUILD_DIR)/bin/NoSpherA2_arm64" \
	  "$(BUILD_DIR)/bin/NoSpherA2_x86_64"

NoSpherA2_arm64:
	@$(MAKE) ARCH=arm64 NoSpherA2
	@cp -f "$(BUILD_DIR)/bin/NoSpherA2" "$(BUILD_DIR)/bin/NoSpherA2_arm64"

NoSpherA2_x86_64:
	@$(MAKE) ARCH=x86_64 NoSpherA2
	@cp -f "$(BUILD_DIR)/bin/NoSpherA2" "$(BUILD_DIR)/bin/NoSpherA2_x86_64"

clean:
	@$(MAKE) -C $(STATIC_DIR) clean
	@$(MAKE) -C $(EXE_DIR) clean
	@$(MAKE) -C $(TEST_DIR) clean
	-@echo ' '

.PHONY: all static tests clean NoSpherA2 lipo NoSpherA2_arm64 NoSpherA2_x86_64
