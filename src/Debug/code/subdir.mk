################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../code/Config.cpp \
../code/Data.cpp \
../code/FreeParam.cpp \
../code/Model.cpp \
../code/Parser.cpp \
../code/Particle.cpp \
../code/Pheromones.cpp \
../code/Swarm.cpp \
../code/Utils.cpp 

OBJS += \
./code/Config.o \
./code/Data.o \
./code/FreeParam.o \
./code/Model.o \
./code/Parser.o \
./code/Particle.o \
./code/Pheromones.o \
./code/Swarm.o \
./code/Utils.o 

CPP_DEPS += \
./code/Config.d \
./code/Data.d \
./code/FreeParam.d \
./code/Model.d \
./code/Parser.d \
./code/Particle.d \
./code/Pheromones.d \
./code/Swarm.d \
./code/Utils.d 


# Each subdirectory must supply rules for building sources it contributes
code/%.o: ../code/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


