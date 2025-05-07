###################################################################################################
# vcfgl Makefile
####################################################################################################

default: all

PROGRAM = vcfgl

# Compiler specs
CXX ?= g++


# No-compile targets
NO_COMPILE = clean test help


VAL_ADD_CRYPTOLIB = -lcrypto
VAL_NOTADD_CRYPTOLIB =

PRINTF_BOLD := $(shell command -v tput > /dev/null && tput bold || echo "")
PRINTF_GREEN := $(shell command -v tput > /dev/null && tput setaf 2 || echo "")
PRINTF_RED := $(shell command -v tput > /dev/null && tput setaf 1 || echo "")
PRINTF_YELLOW := $(shell command -v tput > /dev/null && tput setaf 3 || echo "")
PRINTF_NORMAL := $(shell command -v tput > /dev/null && tput sgr0 || echo "")

####################################################################################################
# [BLOCK START]
# ony run if make will try compiling, i.e. not NO_COMPILE

ifeq (,$(filter $(NO_COMPILE),$(MAKECMDGOALS))) #1

$(info ________________________________________________________________________________)
$(info )
$(info [INFO]    Checking for library availability)
$(info )


## -> [cryptolib availability check]
CRYPTO_TRY = $(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto 2>/dev/null -o /dev/null; echo $$?)

ifeq "$(CRYPTO_TRY)" "0" #1_1

$(info $(PRINTF_GREEN)[INFO]    -> Crypto library is available to link$(PRINTF_NORMAL))
THIS_CRYPTOLIB = $(VAL_ADD_CRYPTOLIB)

else  #1_1

# CRYPTO_TRY != 0

THIS_CRYPTOLIB = $(VAL_NOTADD_CRYPTOLIB)

$(info [INFO]    -> Crypto library is not available to link)
endif #1_1

## -> [htslib source check]

#if htslib source is defined
ifdef HTSSRC #1_2

#if hts source is set to systemwide
ifeq ($(HTSSRC),systemwide) #1_2_1

$(info [INFO]    -> HTSSRC set to systemwide; assuming systemwide installation)
THIS_LIBHTS := -lhts

else #1_2_1

#if hts source path is given
# Adjust $(HTSSRC) to point to your top-level htslib directory
$(info [INFO]    -> HTSSRC is defined as $(HTSSRC))
CPPFLAGS := -I$(realpath $(HTSSRC))
THIS_LIBHTS := $(realpath $(HTSSRC))/libhts.a

endif #1_2_1

#if htssrc not defined
else #1_2

$(info [INFO]    -> HTSSRC is not defined; using htslib submodule)

HTSSRC := $(realpath $(CURDIR)/htslib)
CPPFLAGS := -I$(HTSSRC)
THIS_LIBHTS := $(HTSSRC)/libhts.a

all: .activate_module

endif #1_2

DEV_FLAGS = -g -Wall 
OPTIM_OFF = -O0
OPTIM_ON = -O3

PREV_BUILD_MODE := $(shell grep -oP 'define DEV \K\d' dev.h)


ifeq (dev,$(filter dev,$(MAKECMDGOALS))) #1_3

THIS_BUILD_MODE := 1

ifeq (1,$(PREV_BUILD_MODE)) #1_3_1

# PREV_BUILD_MODE is 1
BUILD_MODE_CHANGED := 0

else #1_3_1

# PREV_BUILD_MODE is 0
BUILD_MODE_CHANGED := 1

endif #1_3_1
else #1_3

THIS_BUILD_MODE := 0

ifeq (1,$(PREV_BUILD_MODE)) #1_3_2
# PREV_BUILD_MODE is 1

BUILD_MODE_CHANGED := 1

else #1_3_2
# PREV_BUILD_MODE is 0

BUILD_MODE_CHANGED := 0

endif #1_3_2
endif #1_3




ifeq (1,$(BUILD_MODE_CHANGED)) #1_4
$(info [INFO]    -> Build mode has been changed)
$(info [INFO]    -> Updating dev.h)
$(shell sed -i 's/define DEV $(PREV_BUILD_MODE)/define DEV $(THIS_BUILD_MODE)/' dev.h)
endif #1_4



ifeq (dev,$(filter dev,$(MAKECMDGOALS))) #1_5

$(info )
$(info ________________________________________________________________________________)
$(info )
$(info [INFO]    Compiling in developer mode)
$(info )


OPTIM_FLAGS := $(OPTIM_OFF)
THIS_MODE_FLAGS := $(DEV_FLAGS) $(OPTIM_FLAGS)


else #1_5

$(info )
$(info ________________________________________________________________________________)
$(info )
$(info [INFO]    Compiling in release mode)
$(info )

OPTIM_FLAGS := $(OPTIM_ON)
THIS_MODE_FLAGS := $(OPTIM_FLAGS)

endif #1_5


$(info [INFO]    -> CXXFLAGS was "$(CXXFLAGS)")
CXXFLAGS += $(THIS_MODE_FLAGS)
$(info $(PRINTF_YELLOW)[INFO]    -> Updated CXXFLAGS to "$(CXXFLAGS)"$(PRINTF_NORMAL))


$(info [INFO]    -> LIBS was "$(LIBS)")
LIBS := $(THIS_LIBHTS) $(THIS_CRYPTOLIB) -lz -lm -lbz2 -llzma -lcurl -lpthread 
$(info $(PRINTF_YELLOW)[INFO]    -> Updated LIBS to "$(LIBS)"$(PRINTF_NORMAL))

$(info )
$(info ________________________________________________________________________________)
$(info )

endif #1

# [BLOCK END]
####################################################################################################


####################################################################################################


dev: all

all: $(PROGRAM)

CXXSRC := $(wildcard *.cpp)

# Preprocessed C++ files
PREP := $(CXXSRC:.cpp=.ii)

# Assembly source files
ASM := $(CXXSRC:.cpp=.s)

# Object files
OBJ := $(CXXSRC:.cpp=.o)

# Dependency files
DEP := $(OBJ:.o=.d)

# don’t pull in .d files for no‐compile targets
ifeq ($(filter $(NO_COMPILE),$(MAKECMDGOALS)),)
-include $(DEP)
endif


####################################################################################################
## [install]
.PHONY: install

prefix          = /usr/local
exec_prefix     = $(prefix)
bindir          = $(exec_prefix)/bin

INSTALL         = install
INSTALL_DIR     = $(INSTALL) -dm0755
INSTALL_PROGRAM = $(INSTALL) -m0755

install: $(PROGRAM)
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(PROGRAM) $(DESTDIR)$(bindir)

####################################################################################################

FLAGS := $(CPPFLAGS) $(CXXFLAGS) 


####################################################################################################
# Versioning
VERSIONNO = 1.2.0

VERSION := v$(VERSIONNO)-$(shell git describe --always --exclude '*')

VERSIONH = version.h

$(VERSIONH):
	@echo '#define VCFGL_VERSION "$(VERSION)"' > $@

####################################################################################################
# Build info
BUILDH = build.h

$(BUILDH): 
	$(info [INFO]    Writing build info to build.h)
	$(shell echo '#define VCFGL_MAKE_CXX ("$(CXX)")' > build.h)
	$(shell echo '#define VCFGL_MAKE_LIBS ("$(LIBS)")' >> build.h)
	$(shell echo '#define VCFGL_MAKE_FLAGS ("$(FLAGS)")' >> build.h)
	$(shell echo '#define VCFGL_MAKE_HTSSRC ("$(HTSSRC)")' >> build.h)
	$(shell echo '#define VCFGL_MAKE_CXXFLAGS ("$(CXXFLAGS)")' >> build.h)
	$(shell echo '#define VCFGL_MAKE_CPPFLAGS ("$(CPPFLAGS)")' >> build.h)


$(PROGRAM): $(OBJ) 
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "$(PRINTF_YELLOW)[INFO]    -> Finishing up$(PRINTF_NORMAL)"
	$(CXX) -o $@ $^ $(LIBS) 
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "$(PRINTF_BOLD)$(PRINTF_GREEN)[FINISHED]    $(PROGRAM) is now ready to use!"
	@echo "              Full path to the program: $(CURDIR)/$(PROGRAM)"
	@echo ""
	@echo "              To get started, run:"
	@echo "              $(CURDIR)/$(PROGRAM) -h"
	@echo "              or:"
	@echo "              ./$(PROGRAM) -h$(PRINTF_NORMAL)"
	@echo ""

%.o: %.cpp $(VERSIONH) $(BUILDH) 
	@echo ""
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "$(PRINTF_YELLOW)[INFO]    -> Compiling $*.cpp$(PRINTF_NORMAL)"
	$(CXX) -c  $(FLAGS) $*.cpp
	$(CXX) -MM $(FLAGS) $*.cpp >$*.d



####################################################################################################
## [clean]
# - clean up the directory

.PHONY: clean
clean:
	$(RM) $(OBJ) $(DEP) $(PREP) $(ASM) $(VERSIONH) $(BUILDH) $(PROGRAM)

####################################################################################################
## [test]
# - run unit tests

.PHONY: test test-%
test:
	bash test/runTests.sh

test-%:
	bash test/runTests.sh $*

####################################################################################################
## [.activate_module]

.PHONY: .activate_module

.activate_module:
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "$(PRINTF_YELLOW)[INFO]	Activating HTSlib submodule$(PRINTF_NORMAL)"
	@echo ""
	git submodule update --init --recursive
	$(MAKE) -C $(HTSSRC)
	@echo "$(PRINTF_GREEN)[INFO]	-> HTSlib submodule is now activated$(PRINTF_NORMAL)"
	@echo ""
	@echo "________________________________________________________________________________"
	

####################################################################################################
## [help]

.PHONY: help
help:
	@echo ""
	@echo "----------------------------------------"
	@echo " Program: $(PROGRAM)"
	@echo " Version: $(VERSION)"
	@echo " License: GNU GPLv3.0"
	@echo "----------------------------------------"
	@echo ""
	@echo " Usage:"
	@echo "   make [target] [FLAG=value...]"
	@echo ""
	@echo " Targets:"
	@echo "   help    - Print this help message"
	@echo "   install - Install the program to $(bindir)"
	@echo "   dev     - Compile in developer/debug mode (activates flags: -g -Wall -O0)"
	@echo "   clean   - Clean up the directory"
	@echo "   test    - Run unit tests"
	@echo ""
	@echo " Flags:"
	@echo "   HTSSRC  - Specifies the source of HTSlib."
	@echo "       Values:"
	@echo "       (empty)          - Use the HTSlib submodule [default]"
	@echo "       systemwide       - Use the systemwide HTSlib installation"
	@echo "       /path/to/htslib  - Use the HTSlib installation at /path/to/htslib"
	@echo ""
	@echo " Examples:"
	@echo "   make                         - Compile in release mode using HTSlib submodule"
	@echo "   make HTSSRC=/path/to/htslib  - Compile in release mode using /path/to/htslib"
	@echo "   make dev HTSSRC=systemwide   - Compile in developer mode using the systemwide HTSlib installation"
	@echo ""
	@echo " Note: If no values are provided for HTSSRC, CXX, CXXFLAGS, or LIBS, defaults will be used."
	@echo ""

####################################################################################################
