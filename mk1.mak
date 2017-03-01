###################################################################
##   THIS MAKE FILE IS GENERATED ONLY FOR TEST AND CUSTOMIZED COMPILE
##   IT COMPILES ONLY A FEW FILES NOT ALL [LATER REMOVEFROM THE PACKAGE]
###################################################################
# Add inputs and outputs from these tool invocations to the build variables 
BINDIR= ./bin
CC=mpic++
CompilerVer=gnu++0x
LogFile=log1.txt



CPP_DEPS = ./bin/src/GenFit2.d \
./bin/src/code/Config.d \
./bin/src/code/Data.d \
./bin/src/code/FreeParam.d \
./bin/src/code/Model.d \
./bin/src/code/Parser.d \
./bin/src/code/Particle.d \
./bin/src/code/Pheromones.d \
./bin/src/code/Swarm.d \
./bin/src/code/Utils.d 

################################################################################
# From File Objects
################################################################################
USER_OBJS := ./lib/libboost_iostreams.a ./lib/libboost_regex.a ./lib/libboost_mpi.a ./lib/libboost_program_options.a ./lib/libboost_filesystem.a ./lib/libboost_system.a ./lib/libboost_serialization.a

LIBS := -lrt

OBJS = ./bin/src/GenFit2.o \
./bin/src/code/Config.o \
./bin/src/code/Data.o \
./bin/src/code/FreeParam.o \
./bin/src/code/Model.o \
./bin/src/code/Parser.o \
./bin/src/code/Particle.o \
./bin/src/code/Pheromones.o \
./bin/src/code/Swarm.o \
./bin/src/code/Utils.o 




################################################################################
# From file: Sources
################################################################################

O_SRCS := 
CPP_SRCS := 
C_UPPER_SRCS := 
C_SRCS := 
S_UPPER_SRCS := 
OBJ_SRCS := 
ASM_SRCS := 
CXX_SRCS := 
C++_SRCS := 
CC_SRCS := 
OBJS := 
C++_DEPS := 
C_DEPS := 
CC_DEPS := 
CPP_DEPS := 
EXECUTABLES := 
CXX_DEPS := 
C_UPPER_DEPS := 

# Every subdirectory with source files must be described here
SUBDIRS := \
src/code \
src \

-include makefile.init
RM := rm -rf




#install:
#	cd $(BINDIR); 



	
#$(CPP_DEPS)   later add
bin/src/GenFit2.o: ./src/GenFit2.cpp ./src/GenFit2.hh 
	pwd
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CC) -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L  -I./include -I/usr/lib/openmpi/include -O3 -g3 -Wall -c -fmessage-length=0 -std=$(CompilerVer) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


#$(CPP_DEPS)   later add
bin/src/code/Utils.o: ./src/code/Utils.cpp ./src/code/Utils.hh 
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CC) -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L  -I./include -I/usr/lib/openmpi/include -O3 -g3 -Wall -c -fmessage-length=0 -std=$(CompilerVer) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
	
	
bin/BioNetFit: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	mpic++ -L/usr/lib -o "BioNetFit" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	
cleanup:
	rm -f ./bin/*.exe
	rm -f ./bin/src/GenFit2.d
	rm -f ./bin/src/GenFit2.o
#	rm -f ./bin/src/code/Utils.o
#	rm -f ./bin/src/code/Utils.d
	