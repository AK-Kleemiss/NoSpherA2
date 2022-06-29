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

all:
	@echo "I think this is $(NAME)"
ifeq ($(NAME),WINDOWS)
	cd Windows && msbuild NoSpherA2.sln
endif
ifeq ($(NAME),LINUX)
	@echo "Start making Linux executable"
	rm -f NoSpherA2
	cd Linux && rm -f NoSpherA2 && make all
endif
ifeq ($(NAME),MAC)
	@echo "Start making Mac executable"
	cd Mac && rm -f NoSpherA2 && make all
endif

test: 
	cd tests && make all -k -B
