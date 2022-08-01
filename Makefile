CXX ?= g++

#  -g     add debugging info to the executable 
#  -Wall  turn on most compiler warnings
CXXFLAGS  = -g -Wall
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
CXXFLAGS := -I"$(HTSSRC)"
LIBHTS := $(HTSSRC)/libhts.a


all: .activate_module

endif

.PHONY: .activate_module test

.activate_module:
	git submodule update --init --recursive
	$(MAKE) -C $(HTSSRC)


VERSION = v0.1


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
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d



$(PROGRAM): version.h $(OBJ) 
	$(CXX) -o $(PROGRAM) *.o $(LIBHTS) $(LIBS) 

clean:
	$(RM) *.o *.d $(PROGRAM) version.h

test: 
	./vcfgl -in test/t1.vcf -out test/t1_pos00_explode0_test -O v -seed 42 -depth 1 -err 0.01 -pos0 0 -explode 0;
	bash -c "diff -I '^##fileDate' -I '^##source vcfgl version:' test/t1_pos00_explode0_test.vcf test/t1_pos00_explode0_god.vcf";
	./vcfgl -in test/t1.vcf -out test/t1_pos00_explode1_test -O v -seed 42 -depth 1 -err 0.01 -pos0 0 -explode 1;
	bash -c "diff -I '^##fileDate' -I '^## source vcfgl version:' test/t1_pos00_explode1_test.vcf test/t1_pos00_explode1_god.vcf";
	./vcfgl -in test/t1.vcf -out test/t1_pos01_explode0_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 0;
	bash -c "diff -I '^##fileDate' -I '^## source vcfgl version:' test/t1_pos01_explode0_test.vcf test/t1_pos01_explode0_god.vcf";
	./vcfgl -in test/t1.vcf -out test/t1_pos01_explode1_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 1;
	bash -c "diff -I '^##fileDate' -I '^## source vcfgl version:' test/t1_pos01_explode1_test.vcf test/t1_pos01_explode1_god.vcf";
	./vcfgl -in test/t2.vcf -out test/t2_pos01_explode0_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 0;
	bash -c "diff -I '^##fileDate' -I '^## source vcfgl version:' test/t2_pos01_explode0_test.vcf test/t2_pos01_explode0_god.vcf";
	./vcfgl -in test/t2.vcf -out test/t2_pos01_explode1_test -O v -seed 42 -depth 1 -err 0.01 -pos0 1 -explode 1;
	bash -c "diff -I '^##fileDate' -I '^## source vcfgl version:' test/t2_pos01_explode1_test.vcf test/t2_pos01_explode1_god.vcf";

