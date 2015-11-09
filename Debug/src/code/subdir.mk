################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/code/Config.cpp \
../src/code/Data.cpp \
../src/code/Model.cpp \
../src/code/Parser.cpp \
../src/code/Particle.cpp \
../src/code/Pheromones.cpp \
../src/code/Swarm.cpp \
../src/code/Utils.cpp \
../src/code/log.cpp 

OBJS += \
./src/code/Config.o \
./src/code/Data.o \
./src/code/Model.o \
./src/code/Parser.o \
./src/code/Particle.o \
./src/code/Pheromones.o \
./src/code/Swarm.o \
./src/code/Utils.o \
./src/code/log.o 

CPP_DEPS += \
./src/code/Config.d \
./src/code/Data.d \
./src/code/Model.d \
./src/code/Parser.d \
./src/code/Particle.d \
./src/code/Pheromones.d \
./src/code/Swarm.d \
./src/code/Utils.d \
./src/code/log.d 


# Each subdirectory must supply rules for building sources it contributes
src/code/%.o: ../src/code/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -I/usr/include/boost -I/usr/include/c++/4.9.2/ -I/usr/lib/openmpi/include -O0 -g3 -Wall -c -fmessage-length=0 -std=c++14 -Wl,-Bstatic -lboost_program_options -Wl,-Bdynamic -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


