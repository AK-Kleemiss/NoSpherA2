REPO_ROOT ?= ../..
CONFIG ?= Release

BUILD_DIR := $(REPO_ROOT)/build/$(CONFIG)

LIBTOOL ?= libtool

SRC_DIR := $(REPO_ROOT)/Src/core
OUT_DIR := $(BUILD_DIR)/obj/core_$(ARCH)

include ../common.mk

# Core sources
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
SRCS := $(filter-out $(SRC_DIR)/mat_nos_math.cpp $(SRC_DIR)/vec_nos_math.cpp, $(SRCS))

OBJ_NAMES := $(basename $(notdir $(SRCS)))
OBJS := $(addprefix $(OUT_DIR)/,$(OBJ_NAMES:=.o))

# Optional PCH (Src/core/pch.h might not exist in all configs)
PCH_HEADER := $(SRC_DIR)/pch.h
PCH_ACTIVE := $(OUT_DIR)/pch.pch
HAVE_PCH := $(wildcard $(PCH_HEADER))

ifeq ($(strip $(HAVE_PCH)),)
  CXXFLAGS := $(CXXFLAGS_COMMON)
else
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

# ---- BasisSetConverter hook ----
$(REPO_ROOT)/Src/basis_data.cpp:
	@echo "Running BasisSetConverter..."
	@$(REPO_ROOT)/Lib/BasisSetConverter $(REPO_ROOT)/BasisSetGenerator/basis_sets $(REPO_ROOT)/Src
	@echo "basis_data.cpp generated"

$(REPO_ROOT)/Src/basis_data.o: $(REPO_ROOT)/Src/basis_data.cpp $(PCH_ACTIVE)
	@echo 'Building file: $<'
	@$(COMP) $(CXXFLAGS) -c "$<" -o "$@"

ifeq ($(strip $(HAVE_PCH)),)
$(REPO_ROOT)/Src/basis_data.o: $(REPO_ROOT)/Src/basis_data.cpp
	@echo 'Building file: $<'
	@$(COMP) $(CXXFLAGS) -c "$<" -o "$@"
endif

OBJS += $(REPO_ROOT)/Src/basis_data.o

LIB_OUT ?= $(BUILD_DIR)/lib/libnosphera2_$(ARCH).a

all: $(LIB_OUT)

$(LIB_OUT): $(OBJS)
	@echo 'Archiving: $@'
	@mkdir -p "$(dir $(LIB_OUT))"
	@$(LIBTOOL) -static -o "$@" $(OBJS)
	@echo 'Done building: $@'

clean:
	-${RM} $(LIB_OUT) $(OBJS) $(OBJS:.o=.d) $(PCH_ACTIVE)
	-@echo ' '

.PHONY: all clean
.SECONDARY:
