OBJS = picoDST.o StRefMultCorr.o  analysis.o 
EXE = analysis



ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++
FLAGS = -m32 -fno-inline -Wall -g  $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX)  -m32 -fno-inline -g -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)

%.o: %.cxx
	$(COMPILE)  $< 

clean:
	rm -f *.o ${EXE} *.root temp_gccflags.c
