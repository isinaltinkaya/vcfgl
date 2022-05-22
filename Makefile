CXX ?= g++

#  -g     add debugging info to the executable 
#  -Wall  turn on most compiler warnings
CXXFLAGS  = -g -Wall
LIBS = -lz -lhts

# HTSSRC := $(CURDIR)/htslib
HTSSRC := /maps/projects/lundbeck/scratch/pfs488/AMOVA/vcfToGlf/htslib


LIBHTS := $(HTSSRC)/libhts.a
LIBS := $(LIBHTS) $(LIBS)

CPPFLAGS += -I$(HTSSRC)



PROGRAM = vcfReader

CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)


-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d



all: $(PROGRAM) $(OBJ)

$(PROGRAM): $(OBJ)
	$(CXX) $(CPPFLAGS) $(LIBS) -o $(PROGRAM) *.o 

clean:
	$(RM) *.o *.d $(PROGRAM)

