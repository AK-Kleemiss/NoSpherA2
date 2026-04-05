ifeq ($(OS),Windows_NT)
    NOS_EXE ?= ../Windows/x64/Release/NoSpherA2.exe
else
    NOS_EXE ?= ../NoSpherA2
endif
export NOS_EXE


PYTEST ?= uv run pytest

PYTEST_ARGS ?= -v

.PHONY: all all_tests

all:
	$(PYTEST) run_test.py $(PYTEST_ARGS)

all_tests: export RUN_FULL_TEST = 1
all_tests:
	$(PYTEST) run_test.py $(PYTEST_ARGS)
