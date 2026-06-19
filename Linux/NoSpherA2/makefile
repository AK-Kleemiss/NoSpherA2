# Thin wrapper: build NoSpherA2 executable on Linux using shared rules
# Compiles only Src/exe and links against the core static library.

# Defaults (can be overridden by caller)
REPO_ROOT ?= ../..
CONFIG ?= Release

# Ensure core sources are *not* pulled in from common.mk defaults.
SRC_DIR := $(REPO_ROOT)/Src/exe

# Reuse the core PCH for the executable, but make sure it's discoverable:
# common.mk adds $(SRC_DIR) to include paths, so also add Src/core.
PCH_HEADER := $(REPO_ROOT)/Src/core/pch.h
USER_INC_DIRS := $(REPO_ROOT)/Src/core $(REPO_ROOT)/Src/exe $(REPO_ROOT)/Src $(REPO_ROOT)/Lib $(REPO_ROOT)/Src/occ

INCLUDE_BASIS_DATA := 0
SRCS := $(wildcard $(SRC_DIR)/*.cpp)

OUT_DIR = $(BUILD_DIR)/obj/exe

include ../common.mk

OBJ_NAMES := $(basename $(notdir $(SRCS)))
OBJS := $(addprefix $(OUT_DIR)/,$(OBJ_NAMES:=.o))

PCH_GCH_ACTIVE := $(OUT_DIR)/pch.gch
PCH_BASENAME := $(notdir $(PCH_HEADER))
GCC_OPTS_NO_PCH := $(filter-out -include $(PCH_BASENAME),$(GCC_OPTS))

$(PCH_GCH_ACTIVE): $(PCH_HEADER)
	@mkdir -p "$(OUT_DIR)"
	@echo 'Building precompiled header...'
	@$(WITH_ONEAPI) ${COMP} $(GCC_OPTS_NO_PCH) -I"$$MKLROOT/include" $(INCLUDES) -x c++-header -o "$@" "$<"

$(OUT_DIR)/%.o: $(SRC_DIR)/%.cpp $(PCH_GCH_ACTIVE)
	@echo 'Building file: $<'
	@mkdir -p "$(OUT_DIR)"
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS} -I"$$MKLROOT/include" $(INCLUDES) \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

-include $(OBJS:.o=.d)

# Link inputs
LIB_NOS ?= $(BUILD_DIR)/lib/libnosphera2.a

LIBS := -Wl,-Bstatic -lstdc++ -Wl,-Bdynamic -L$(REPO_ROOT)/Lib/occ/lib -ltbb -Wl,-rpath,'$$ORIGIN' -lpthread -lm -lc -lquadmath

ifeq ($(CONFIG),Debug)
all: $(BUILD_DIR)/bin/NoSpherA2_Debug
else
all: $(BUILD_DIR)/bin/NoSpherA2
endif

NoSpherA2: $(BUILD_DIR)/bin/NoSpherA2

$(BUILD_DIR)/bin/NoSpherA2: $(OBJS) $(LIB_NOS) $(PCH_GCH_ACTIVE)
	@echo 'Building target: $@'
	@mkdir -p "$(BUILD_DIR)/bin"
	@$(WITH_ONEAPI) ${COMP} -o "$@" $(OBJS) $(LIB_NOS) ${LIBRASCALINE} $(MKL_LINK_FLAGS) $(LIBCINT) $(REPO_ROOT)/Lib/occ/lib/libocc.a ${LIBS}
	@echo 'Done building NoSpherA2'

NoSpherA2_Debug: $(BUILD_DIR)/bin/NoSpherA2_Debug

$(BUILD_DIR)/bin/NoSpherA2_Debug: $(OBJS) $(LIB_NOS) $(PCH_GCH_ACTIVE)
	@echo 'Building target: $@'
	@mkdir -p "$(BUILD_DIR)/bin"
	@$(WITH_ONEAPI) ${COMP} -o "$@" $(OBJS) $(LIB_NOS) ${LIBRASCALINE} $(MKL_LINK_FLAGS) $(LIBCINT) $(REPO_ROOT)/Lib/occ/lib/libocc.a ${LIBS}
	@echo 'Done building Debug version of NoSpherA2'

clean:
	-${RM} $(BUILD_DIR)/bin/NoSpherA2 $(BUILD_DIR)/bin/NoSpherA2_Debug $(OBJS) $(OBJS:.o=.d) $(PCH_GCH_ACTIVE)
	-@echo ' '

.PHONY: all clean
.SECONDARY:
