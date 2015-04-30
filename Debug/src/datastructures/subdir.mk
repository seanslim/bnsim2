################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/datastructures/common.cpp \
../src/datastructures/moleculeInfo.cpp \
../src/datastructures/spacegrid.cpp 

OBJS += \
./src/datastructures/common.o \
./src/datastructures/moleculeInfo.o \
./src/datastructures/spacegrid.o 

CPP_DEPS += \
./src/datastructures/common.d \
./src/datastructures/moleculeInfo.d \
./src/datastructures/spacegrid.d 


# Each subdirectory must supply rules for building sources it contributes
src/datastructures/%.o: ../src/datastructures/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


