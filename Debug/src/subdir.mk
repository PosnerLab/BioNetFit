################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/GenFit2.cpp 

OBJS += \
./src/GenFit2.o 

CPP_DEPS += \
./src/GenFit2.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -I/usr/include/boost -I/usr/include/c++/4.9.2/ -I/usr/lib/openmpi/include -O0 -g3 -Wall -c -fmessage-length=0 -std=c++14 -Wl,-Bstatic -lboost_program_options -Wl,-Bdynamic -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


