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


SRC = $(wildcard *.c)
CXXSRC = $(wildcard *.cpp)

OBJ := $(CXXSRC:.cpp=.o) $(SRC:.c=.o)


PROGRAMS = $(SRC:.c=) $(CXXSRC:.cpp=)
all: $(PROGRAMS)


FLAGS := $(FLAGS) $(CXXFLAGS)

%: %.cpp
	@echo
	@echo Building $@
	@echo
	$(CXX) -c -g -Wall $(FLAGS) $*.cpp
	$(CXX) -MM  -g -Wall $(FLAGS) $*.cpp > $*.d
	$(CXX) -o $* $*.o $(LIBHTS) -lz -lm -lbz2 -llzma -lcurl -lpthread

%: %.c
	@echo
	@echo Building $@
	@echo
	$(CXX) -c -g -Wall $(FLAGS) $*.c
	$(CXX) -MM  -g -Wall $(FLAGS) $*.c > $*.d
	$(CXX) -o $* $*.o $(LIBHTS) -lz -lm -lbz2 -llzma -lcurl -lpthread


clean:
	$(RM) *.o *.d $(PROGRAMS) 

test:
	$(RM) -r testwd/;
	mkdir -pv testwd;
	$(info Test 1)
	./gtDiscordance -i data/call.vcf -t data/truth.vcf -o testwd/out_test.tsv;
	diff -s testwd/out_test.tsv reference/test_doGq0.tsv
	./gtDiscordance -i data/call3.vcf -t data/truth3.vcf -o testwd/out_test3.tsv -doGQ 7;
	diff -s testwd/out_test3.tsv reference/test3_doGq7.tsv
	./gtDiscordance -i data/call3_nogq_withpl.vcf -t data/truth3.vcf -o testwd/out_test3_nogq_withpl.tsv -doGQ 8;
	diff -s testwd/out_test3_nogq_withpl.tsv reference/test3_doGq7.tsv
	./gtDiscordance -i data/call3_nogq_withpl_unobservedAllele.vcf -t data/truth3.vcf -o testwd/out_test3_nogq_withpl_unobservedAllele.tsv -doGQ 8;
	diff -s testwd/out_test3_nogq_withpl_unobservedAllele.tsv reference/test3_doGq7.tsv
	./fetchGl -i data/test12.vcf -gt CC > testwd/test12.csv 2>/dev/null
	diff -s testwd/test12.csv reference/test12.csv
