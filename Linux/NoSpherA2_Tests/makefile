# Thin wrapper: build NoSpherA2 tests on Linux.
# Compiles only Src/tests and links against the core static library.
# Note: The existing test sources use the Visual Studio CppUnitTest framework;
# you'll likely adapt these to a Linux-friendly framework later.

REPO_ROOT ?= ../..

SRC_DIR := $(REPO_ROOT)/Src/tests
SRCS := $(wildcard $(SRC_DIR)/*.cpp)

OUT_DIR = $(BUILD_DIR)/obj/tests

include ../common.mk

# Catch2 is installed via apt on Linux. The Debian/Ubuntu package installs headers
# under /usr/include/catch2.
USER_INC_DIRS += /usr/include

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

$(OUT_DIR)/%.o: $(SRC_DIR)/%.cpp $(PCH_GCH_ACTIVE)
	@echo 'Building file: $<'
	@mkdir -p "$(OUT_DIR)"
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS} -I"$$MKLROOT/include" $(INCLUDES) \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

-include $(OBJS:.o=.d)

LIB_NOS ?= $(BUILD_DIR)/lib/libnosphera2.a

LIBS := -Wl,-Bstatic -lstdc++ -Wl,-Bdynamic -L$(REPO_ROOT)/Lib/occ/lib -ltbb -Wl,-rpath,'$$ORIGIN' -lpthread -lm -lc -lquadmath

all: NoSpherA2_Tests

NoSpherA2_Tests: $(OBJS) $(LIB_NOS) $(PCH_GCH_ACTIVE)
	@echo 'Building target: $@'
	@mkdir -p "$(BUILD_DIR)/bin"
	@$(WITH_ONEAPI) ${COMP} -o "$(BUILD_DIR)/bin/NoSpherA2_Tests" $(OBJS) $(LIB_NOS) ${LIBRASCALINE} $(MKL_LINK_FLAGS) $(LIBCINT) $(REPO_ROOT)/Lib/occ/lib/libocc.a ${LIBS}
	@echo 'Done building: $@'

clean:
	-${RM} $(BUILD_DIR)/bin/NoSpherA2_Tests $(OBJS) $(OBJS:.o=.d) $(PCH_GCH_ACTIVE)
	-@echo ' '

.PHONY: all clean
.SECONDARY:
