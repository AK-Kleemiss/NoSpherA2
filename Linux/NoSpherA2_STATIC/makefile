REPO_ROOT ?= ../..
CONFIG ?= Release

BUILD_DIR := $(REPO_ROOT)/build/$(CONFIG)

# If we're going to link this static library into a shared object, it must be PIC.
ifeq ($(BUILD_SHARED),1)
PIC ?= 1
endif

AR ?= ar  # GNU ar is used for static library creation; override if using a different tool (e.g. lib on Windows)
RANLIB ?= ranlib # ranlib is used to index static libraries on some platforms; override with `true` if not needed (e.g. Windows .lib)

INCLUDE_BASIS_DATA := 1

SRC_DIR := $(REPO_ROOT)/Src/core

OUT_DIR = $(BUILD_DIR)/obj/core

include ../common.mk

# Core sources
SRCS := $(wildcard $(SRC_CORE_DIR)/*.cpp)
SRCS := $(filter-out $(SRC_CORE_DIR)/mat_nos_math.cpp $(SRC_CORE_DIR)/vec_nos_math.cpp, $(SRCS))

OBJ_NAMES := $(basename $(notdir $(SRCS)))
OBJS := $(addprefix $(OUT_DIR)/,$(OBJ_NAMES:=.o))

PCH_HEADER := $(SRC_CORE_DIR)/pch.h
PCH_BASENAME := $(notdir $(PCH_HEADER))
PCH_GCH_ACTIVE := $(OUT_DIR)/pch.gch
GCC_OPTS_NO_PCH := $(filter-out -include $(PCH_BASENAME),$(GCC_OPTS))

$(PCH_GCH_ACTIVE): $(PCH_HEADER)
	@mkdir -p "$(OUT_DIR)"
	@echo 'Building precompiled header...'
	@$(WITH_ONEAPI) ${COMP} $(GCC_OPTS_NO_PCH) -I"$$MKLROOT/include" $(INCLUDES) -x c++-header -o "$@" "$<"

$(OUT_DIR)/%.o: $(SRC_CORE_DIR)/%.cpp $(PCH_GCH_ACTIVE)
	@echo 'Building file: $<'
	@mkdir -p "$(OUT_DIR)"
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS} -I"$$MKLROOT/include" $(INCLUDES) \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

-include $(OBJS:.o=.d)

# ---- BasisSetConverter hook ----
# Actually run the converter to generate the file
$(SRC_CORE_DIR)/../basis_data.cpp:
	@echo "Running BasisSetConverter..."
	@$(REPO_ROOT)/Lib/BasisSetConverter $(REPO_ROOT)/BasisSetGenerator/basis_sets $(REPO_ROOT)/Src
	@echo "basis_data.cpp generated"

$(SRC_CORE_DIR)/../basis_data.o: $(SRC_CORE_DIR)/../basis_data.cpp $(PCH_GCH_ACTIVE)
	@echo 'Building file: $<'
	@mkdir -p "$(OUT_DIR)"
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS} -I"$$MKLROOT/include" $(INCLUDES) \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

OBJS += $(SRC_CORE_DIR)/../basis_data.o

LIB_OUT ?= $(BUILD_DIR)/lib/libnosphera2.a

all: $(LIB_OUT)

$(LIB_OUT): $(OBJS)
	@echo 'Archiving: $@'
	@mkdir -p "$(dir $(LIB_OUT))"
	@$(AR) rcs "$@" $(OBJS)
	@$(RANLIB) "$@" 2>/dev/null || true
	@echo 'Done building: $@'

clean:
	-${RM} $(LIB_OUT) $(OBJS) $(OBJS:.o=.d) $(PCH_GCH_ACTIVE)
	-@echo ' '

.PHONY: all clean
.SECONDARY: