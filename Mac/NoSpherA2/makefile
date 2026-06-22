REPO_ROOT ?= ../..
CONFIG ?= Release

BUILD_DIR := $(REPO_ROOT)/build/$(CONFIG)

SRC_DIR := $(REPO_ROOT)/Src/exe
SRCS := $(wildcard $(SRC_DIR)/*.cpp)

OUT_DIR := $(BUILD_DIR)/obj/exe_$(ARCH)_$(DEPLOYMENT_TAG)

include ../common.mk

OBJ_NAMES := $(basename $(notdir $(SRCS)))
OBJS := $(addprefix $(OUT_DIR)/,$(OBJ_NAMES:=.o))

PCH_HEADER := $(REPO_ROOT)/Src/core/pch.h
PCH_ACTIVE := $(OUT_DIR)/pch.pch
HAVE_PCH := $(wildcard $(PCH_HEADER))

ifeq ($(strip $(HAVE_PCH)),)
  CXXFLAGS := $(CXXFLAGS_COMMON)
else
  # ensure Src/core is on include path
  CXXFLAGS_COMMON := -I$(REPO_ROOT)/Src/core $(CXXFLAGS_COMMON)
  CXXFLAGS := $(CXXFLAGS_COMMON) -include-pch $(PCH_ACTIVE)
endif

$(PCH_ACTIVE): $(PCH_HEADER)
	@mkdir -p "$(OUT_DIR)"
	@echo 'Building precompiled header...'
	@$(COMP) $(CXXFLAGS_COMMON) -x c++-header -o "$@" "$<"

$(OUT_DIR)/%.o: $(SRC_DIR)/%.cpp $(PCH_ACTIVE)
	@echo 'Building file: $<'
	@mkdir -p "$(OUT_DIR)"
	@$(COMP) $(CXXFLAGS) -c "$<" -o "$@"

ifeq ($(strip $(HAVE_PCH)),)
$(OUT_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo 'Building file: $<'
	@mkdir -p "$(OUT_DIR)"
	@$(COMP) $(CXXFLAGS) -c "$<" -o "$@"
endif

-include $(OBJS:.o=.d)

LIB_NOS ?= $(BUILD_DIR)/lib/libnosphera2_$(ARCH).a

ifeq ($(CONFIG),Debug)
all: $(BUILD_DIR)/bin/NoSpherA2_Debug
else
all: $(BUILD_DIR)/bin/NoSpherA2
endif

$(BUILD_DIR)/bin/NoSpherA2: $(OBJS) $(LIB_NOS)
	@echo 'Linking target: $@'
	@mkdir -p "$(BUILD_DIR)/bin"
	@$(COMP) -o "$@" $(OBJS) $(LIB_NOS) $(LIBRASCALINE) $(LIBCINT) $(LIBOCC) $(LDFLAGS_COMMON)
	@echo 'Done building NoSpherA2'

$(BUILD_DIR)/bin/NoSpherA2_Debug: $(OBJS) $(LIB_NOS)
	@echo 'Linking target: $@'
	@mkdir -p "$(BUILD_DIR)/bin"
	@$(COMP) -o "$@" $(OBJS) $(LIB_NOS) $(LIBRASCALINE) $(LIBCINT) $(LIBOCC) $(LDFLAGS_COMMON)
	@echo 'Done building NoSpherA2 (Debug)'

clean:
	-${RM} $(BUILD_DIR)/bin/NoSpherA2 $(BUILD_DIR)/bin/NoSpherA2_Debug $(OBJS) $(OBJS:.o=.d) $(PCH_ACTIVE)
	-@echo ' '

.PHONY: all clean
.SECONDARY:
