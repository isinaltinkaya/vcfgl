CXX ?= g++


ifeq ($(DEV),1)
$(info Compiling in developer mode)
COMPILEMODE = -DDEV=1
# TODO disable clean and replace with ngsAMOVA style checks
all: clean
CXXFLAGS := -g -Wall -O0
else
$(info Compiling in release mode; will enable optimizaton (-O3))
COMPILEMODE = -DDEV=0
# TODO disable clean and replace with ngsAMOVA style checks
all: clean
CXXFLAGS := -O3
endif

LIBS = -lz -lm -lbz2 -llzma -lcurl -lpthread


#if htslib source is defined
ifdef HTSSRC

#if hts source is set to systemwide
ifeq ($(HTSSRC),systemwide)
$(info HTSSRC set to systemwide; assuming systemwide installation)
LIBHTS := -lhts

else

#if hts source path is given
# Adjust $(HTSSRC) to point to your top-level htslib directory
$(info HTSSRC defined: $(HTSSRC))
CXXFLAGS += -I"$(realpath $(HTSSRC))"
LIBHTS := $(realpath $(HTSSRC))/libhts.a

endif

#if htssrc not defined
else

$(info HTSSRC not defined; using htslib submodule)
$(info Use `make HTSSRC=/path/to/htslib` to build using a local htslib installation)
$(info Use `make HTSSRC=systemwide` to build using the systemwide htslib installation)

HTSSRC := $(realpath $(CURDIR)/htslib)
CXXFLAGS += -I"$(HTSSRC)"
LIBHTS := $(HTSSRC)/libhts.a


all: .activate_module

endif

.PHONY: .activate_module test 

.activate_module:
	git submodule update --init --recursive
	$(MAKE) -C $(HTSSRC)


VERSION = v0.3.1


ifneq "$(wildcard .git)" ""
	VERSION += $(shell git describe --always --dirty)
endif


version.h:
	echo '#define VCFGL_VERSION "$(VERSION)"' > $@



PROGRAM = vcfgl
all: $(PROGRAM)

CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $(COMPILEMODE) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $(COMPILEMODE) $*.cpp >$*.d



$(PROGRAM): version.h $(OBJ) 
	$(CXX) -o $(PROGRAM) *.o $(LIBHTS) $(LIBS) 

clean:
	$(RM) *.o *.d $(PROGRAM) version.h

test: 
	bash test/runTests.sh
	./vcfgl

