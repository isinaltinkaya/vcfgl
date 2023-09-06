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


VERSION = v0.2


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
	./vcfgl -in test/t1.vcf -out test/t1_pos00_explode0_test -O v -seed 42 -depth 1 -err 0.01 -pos0 0 -explode 0;
	bash -c "diff -I '^##'  test/t1_pos00_explode0_test.vcf test/reference/t1_pos00_explode0.vcf";
	./vcfgl -in test/t1.vcf -out test/t1_pos00_explode1_test -O v -seed 42 -depth 1 -err 0.01 -pos0 0 -explode 1;
	bash -c "diff -I '^##'  test/t1_pos00_explode1_test.vcf test/reference/t1_pos00_explode1.vcf";
	./vcfgl -in test/t1.vcf -out test/t1_pos01_explode0_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 0;
	bash -c "diff -I '^##'  test/t1_pos01_explode0_test.vcf test/reference/t1_pos01_explode0.vcf";
	./vcfgl -in test/t1.vcf -out test/t1_pos01_explode1_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 1;
	bash -c "diff -I '^##'  test/t1_pos01_explode1_test.vcf test/reference/t1_pos01_explode1.vcf";
	./vcfgl -in test/t2.vcf -out test/t2_pos01_explode0_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 0;
	bash -c "diff -I '^##'  test/t2_pos01_explode0_test.vcf test/reference/t2_pos01_explode0.vcf";
	./vcfgl -in test/t2.vcf -out test/t2_pos01_explode1_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 1;
	bash -c "diff -I '^##'  test/t2_pos01_explode1_test.vcf test/reference/t2_pos01_explode1.vcf";
	./vcfgl -in test/t2.vcf -out test/t2_pos01_explode1_gp_gl_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 1 -addGP 1 -addPL 1;
	bash -c "diff -I '^##'  test/t2_pos01_explode1_gp_gl_test.vcf test/reference/t2_pos01_explode1_gp_gl.vcf";
	./vcfgl -in test/t2.vcf -out test/t2_pos01_explode1_gp_gl_qs_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 1 -addQS 1 -addGP 1 -addPL 1;
	bash -c "diff -I '^##'  test/t2_pos01_explode1_gp_gl_qs_test.vcf test/reference/t2_pos01_explode1_gp_gl_qs.vcf";

