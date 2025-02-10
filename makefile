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

NoSpherA2_Debug:
	@echo "Building NoSpherA2_Debug for $(NAME)"
ifeq ($(NAME),WINDOWS)
	cd Windows && msbuild NoSpherA2.sln /p:Configuration=Debug /p:Platform=x64 && cd .. && copy Windows\x64\Debug\NoSpherA2.exe .
endif
ifeq ($(NAME),LINUX)
	@echo "Start making Linux executable for NoSpherA2_Debug"
	rm -f NoSpherA2_Debug
	cd Linux && rm -f NoSpherA2_Debug && make NoSpherA2_Debug -j
endif
ifeq ($(NAME),MAC)
	@echo "Start making Mac executable for NoSpherA2_Debug"
	cd Mac && rm -f NoSpherA2_Debug && make NoSpherA2_Debug -j
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


.PHONY: test tests NoSpherA2 all NoSpherA2_Debug
