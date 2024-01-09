CXX ?= g++

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

$(info HTSSRC not defined; will try to use main program submodule)

HTSSRC := $(realpath $(CURDIR)/../htslib)
CXXFLAGS += -I"$(HTSSRC)"
LIBHTS := $(HTSSRC)/libhts.a


endif

.PHONY: test 


PROGRAM = gtDiscordance


all: $(PROGRAM)

CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

FLAGS := $(FLAGS) $(CXXFLAGS)

-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c  $(FLAGS) $*.cpp
	$(CXX) -MM $(FLAGS) $*.cpp >$*.d



$(PROGRAM): $(OBJ) 
	$(CXX) -o $(PROGRAM) *.o $(LIBHTS) $(LIBS) 

clean:
	$(RM) *.o *.d $(PROGRAM) 

test:
	$(RM) -r testwd/;
	mkdir -pv testwd;
	$(info Test 1)
	./$(PROGRAM) -i data/call.vcf -t data/truth.vcf -o testwd/out_test.tsv;
	diff -s testwd/out_test.tsv reference/test_doGq0.tsv
	./$(PROGRAM) -i data/call3.vcf -t data/truth3.vcf -o testwd/out_test3.tsv -doGQ 7;
	diff -s testwd/out_test3.tsv reference/test3_doGq7.tsv

