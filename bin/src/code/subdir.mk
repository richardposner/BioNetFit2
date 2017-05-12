################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/code/Config.cpp \
../src/code/Data.cpp \
../src/code/FreeParam.cpp \
../src/code/Model.cpp \
../src/code/Parser.cpp \
../src/code/Particle.cpp \
../src/code/Pheromones.cpp \
../src/code/Evaluate.cpp \
../src/code/Swarm.cpp \
../src/code/Utils.cpp 

OBJS += \
./src/code/Config.o \
./src/code/Data.o \
./src/code/FreeParam.o \
./src/code/Model.o \
./src/code/Parser.o \
./src/code/Particle.o \
./src/code/Pheromones.o \
./src/code/Evaluate.o \
./src/code/Swarm.o \
./src/code/Utils.o 

CPP_DEPS += \
./src/code/Config.d \
./src/code/Data.d \
./src/code/FreeParam.d \
./src/code/Model.d \
./src/code/Parser.d \
./src/code/Particle.d \
./src/code/Pheromones.d \
./src/code/Evaluate.d \
./src/code/Swarm.d \
./src/code/Utils.d 


# Each subdirectory must supply rules for building sources it contributes
src/code/%.o: ../src/code/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -w -fpermissive -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -I../include -I/usr/lib/openmpi/include -O3 -g3 -Wall -c -fmessage-length=0 -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


