CC      ?= gcc
CXX     ?= g++

SRC     := $(wildcard *.c)
CXXSRC  := $(wildcard *.cpp)
OBJ     := $(CXXSRC:.cpp=.o) $(SRC:.c=.o)
PROGRAMS := $(SRC:.c=) $(CXXSRC:.cpp=)

LIBS    := -lz -lm -lbz2 -llzma -lcurl -lpthread

ifdef HTSSRC
    ifeq ($(HTSSRC),systemwide)
        $(info HTSSRC set to systemwide; assuming systemwide installation)
        LIBHTS := -lhts
    else
        $(info HTSSRC defined: $(HTSSRC))
        CXXFLAGS += -I"$(realpath $(HTSSRC))"
        LIBHTS := $(realpath $(HTSSRC))/libhts.a
    endif
else
    $(info HTSSRC not defined; will try to use main program submodule)
    HTSSRC := $(realpath $(CURDIR)/../htslib)
    CXXFLAGS += -I"$(HTSSRC)"
    LIBHTS := $(HTSSRC)/libhts.a
endif

# === Build configurations ===

DEV_WARNINGS := -Wall -Wextra -Werror -pedantic -Wshadow -Wconversion -Wundef -Wstrict-prototypes -Wpointer-arith -Wcast-align
DEV_FLAGS    := -std=c17 $(DEV_WARNINGS) -g -O0 -fsanitize=address,undefined
RELEASE_FLAGS := -std=c17 -O3 -DNDEBUG


ifeq ($(MAKECMDGOALS),dev)
    BUILD_FLAGS := $(DEV_FLAGS)
endif

CFLAGS   := $(BUILD_FLAGS)
CXXFLAGS := $(BUILD_FLAGS)
LDFLAGS  := $(filter -fsanitize=%, $(BUILD_FLAGS))  # extract sanitizer flags if any

# ========= Targets =========

all: $(PROGRAMS)

dev: $(PROGRAMS)

%: %.cpp
	@echo
	@echo "Building $@ [$(@:dev=%)]"
	@echo
	$(CXX) -c $(CXXFLAGS) $< -o $*.o
	$(CXX) -MM $(CXXFLAGS) $< > $*.d
	$(CXX) $(LDFLAGS) -o $@ $*.o $(LIBHTS) $(LIBS)

%: %.c
	@echo
	@echo "Building $@ [$(@:dev=%)]"
	@echo
	$(CC) -c $(CFLAGS) $< -o $*.o
	$(CC) -MM $(CFLAGS) $< > $*.d
	$(CC) $(LDFLAGS) -o $@ $*.o $(LIBHTS) $(LIBS)

clean:
	$(RM) *.o *.d $(PROGRAMS)

.PHONY: all clean dev test

test:
	$(RM) -r testwd/
	mkdir -pv testwd
	@echo
	@echo "Running gtDiscordance test 1"
	./gtDiscordance -i data/call.vcf -t data/truth.vcf -o testwd/out_test.tsv
	diff -s testwd/out_test.tsv reference/test_doGq0.tsv
	@echo
	@echo "Running gtDiscordance test 2"
	./gtDiscordance -i data/call3.vcf -t data/truth3.vcf -o testwd/out_test3.tsv -doGQ 7
	diff -s testwd/out_test3.tsv reference/test3_doGq7.tsv
	@echo
	@echo "Running gtDiscordance test 3"
	./gtDiscordance -i data/call3_nogq_withpl.vcf -t data/truth3.vcf -o testwd/out_test3_nogq_withpl.tsv -doGQ 8
	diff -s testwd/out_test3_nogq_withpl.tsv reference/test3_doGq7.tsv
	@echo
	@echo "Running gtDiscordance test 4"
	./gtDiscordance -i data/call3_nogq_withpl_unobservedAllele.vcf -t data/truth3.vcf -o testwd/out_test3_nogq_withpl_unobservedAllele.tsv -doGQ 8
	diff -s testwd/out_test3_nogq_withpl_unobservedAllele.tsv reference/test3_doGq7.tsv
	@echo
	@echo "Running fetchGl test 1"
	./fetchGl -i data/test12.vcf -gt CC > testwd/test12.csv 2>/dev/null
	diff -s testwd/test12.csv reference/test12.csv
	@echo
	@echo "Running setAlleles test 1"
	./setAlleles -i data/data2.vcf -a data/data2_alleles1.tsv -o testwd/data2_alleles1 -O v
	diff -s testwd/data2_alleles1.vcf reference/data2_alleles1.vcf
	@echo
	@echo "Running setAlleles test 2"
	./setAlleles -i data/data2.vcf -a data/data2_alleles2.tsv -o testwd/data2_alleles2 -O v
	diff -s testwd/data2_alleles2.vcf reference/data2_alleles2.vcf
	@echo
	@echo "All tests passed!"
	@echo
	
	
