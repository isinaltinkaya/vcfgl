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



PROGRAM = vcfgl
all: $(PROGRAM)

CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d



$(PROGRAM): $(OBJ)
	$(CXX) -o $(PROGRAM) *.o $(LIBHTS) $(LIBS) 

clean:
	$(RM) *.o *.d $(PROGRAM)

test: 
	./vcfgl -in test/t1.vcf -out test/t1_testvcfgl -err 0.01 -seed 42 -depth 1 -isSim 1 -mode v;
	bash -c "diff -I '^##fileDate' test/t1_vcf_vcfgl.vcf test/t1_testvcfgl.vcf" ;
