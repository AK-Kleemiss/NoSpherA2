ifeq ($(OS), Windows_NT)
  NAME := WINDOWS
else
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S), Linux)
    NAME := LINUX
  endif
  ifeq ($(UNAME_S), Darwin)
    NAME := MAC
  endif
endif

all: NoSpherA2

NoSpherA2_Rascaline:
	@echo "Building NoSpherA2_Rascaline for $(NAME)"
ifeq ($(NAME),WINDOWS)
	cd Windows && msbuild NoSpherA2.sln /p:Configuration=Release /p:Platform=x64 && cd .. && copy Windows\x64\Release\NoSpherA2.exe .
endif
ifeq ($(NAME),LINUX)
	@echo "Start making Linux executable for NoSpherA2_Rascaline"
	rm -f NoSpherA2_Rascaline
	cd Linux && rm -f NoSpherA2_Rascaline && make NoSpherA2_Rascaline -j
endif
ifeq ($(NAME),MAC)
	@echo "Start making Mac executable for NoSpherA2_Rascaline"
	cd Mac && rm -f NoSpherA2_Rascaline && make NoSpherA2_Rascaline -j
endif

NoSpherA2:
	@echo "Building NoSpherA2 for $(NAME)"
ifeq ($(NAME),WINDOWS)
	cd Windows && msbuild NoSpherA2.sln /p:Configuration=Release /p:Platform=x64 && cd .. && copy Windows\x64\Release\NoSpherA2.exe .
endif
ifeq ($(NAME),LINUX)
	@echo "Start making Linux executable"
	rm -f NoSpherA2
	cd Linux && rm -f NoSpherA2 && make all -j
endif
ifeq ($(NAME),MAC)
	@echo "Start making Mac executable"
	cd Mac && rm -f NoSpherA2 && make all -j
endif

test: 
	cd tests && make all -k -B
tests: 
	cd tests && make all -k -B

test_r: 
	@echo "Running tests for Rascaline"
	cd tests && make all -k -B MODE=RASCALINE
tests_r:
	@echo "Running tests for Rascaline"
	cd tests && make all -k -B MODE=RASCALINE


.PHONY: test tests NoSpherA2_Rascaline NoSpherA2 all
