CXX=mpic++
CFLAGS=-Wall -g -I../include -D__GXX_EXPERIMENTAL_CXX0X__ -std=c++0x
LFLAGS=-Wl,-Bstatic -L../lib -lboost_iostreams -lboost_regex -lboost_mpi -lboost_program_options -lboost_filesystem -lboost_system -lboost_serialization -Wl,-Bdynamic -pthread -lrt

.PHONY: all clean

HFILES=\
	GenFit2.hh\
	code/Config.hh\
	code/Evaluate.hh\
	code/Model.hh\
	code/Particle.hh\
	code/Setting.hh\
	code/Timer.hh\
	code/Data.hh\
	code/FreeParam.hh\
	code/Parser.hh\
	code/Pheromones.hh\
	code/Swarm.hh\
	code/Utils.hh

OFILES=\
	code/Config.o\
	code/Data.o\
	code/FreeParam.o\
	code/Model.o\
	code/Parser.o\
	code/Particle.o\
	code/Pheromones.o\
	code/Evaluate.o\
	code/Swarm.o\
	code/Utils.o\
	GenFit2.o

EXES=BioNetFit2

all: $(EXES)
clean:
	rm -rf ../bin/BioNetFit2 $(OFILES) $(EXES)

docker_clean:
	rm $(OFILES)

BioNetFit2: $(HFILES) $(OFILES) Makefile
	$(CXX) -o ../bin/BioNetFit2 $(OFILES) $(LFLAGS)

%.o: %.cpp $(HFILES)
	$(CXX) $(CFLAGS) -o $*.o -c $*.cpp
