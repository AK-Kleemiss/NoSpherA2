NOS_EXE ?= ../NoSpherA2

ifeq ($(OS),Windows_NT)
    PYTEST_CHECK := python -m pytest --version >nul 2>&1
    UV_INSTALL   := powershell -ExecutionPolicy Bypass -Command "irm https://astral.sh/uv/install.ps1 | iex"
    UV_EXE       := $(APPDATA)\uv\bin\uv.exe
else
    PYTEST_CHECK := python3 -m pytest --version >/dev/null 2>&1
    UV_INSTALL   := curl -LsSf https://astral.sh/uv/install.sh | sh
    UV_EXE       := $(HOME)/.local/bin/uv
endif

# Resolve pytest runner once
define PYTEST_RUNNER
$(shell \
    if $(PYTEST_CHECK); then \
        echo python -m pytest; \
    elif command -v uv >/dev/null 2>&1; then \
        echo uv run pytest; \
    elif [ -x "$(UV_EXE)" ]; then \
        echo $(UV_EXE) run pytest; \
    else \
        $(UV_INSTALL) >&2 && echo $(UV_EXE) run pytest; \
    fi \
)
endef

PYTEST ?= $(PYTEST_RUNNER)

PYTEST_ARGS ?= -v

.PHONY: all all_tests

all:
	NOS_EXE=$(NOS_EXE) $(PYTEST) run_test.py $(PYTEST_ARGS)

all_tests:
	NOS_EXE=$(NOS_EXE) RUN_FULL_TEST=1 $(PYTEST) run_test.py $(PYTEST_ARGS)
