################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/agents/ChemotacticBacteria.cpp \
../src/agents/EPS.cpp \
../src/agents/QSBacteria.cpp \
../src/agents/agent.cpp 

OBJS += \
./src/agents/ChemotacticBacteria.o \
./src/agents/EPS.o \
./src/agents/QSBacteria.o \
./src/agents/agent.o 

CPP_DEPS += \
./src/agents/ChemotacticBacteria.d \
./src/agents/EPS.d \
./src/agents/QSBacteria.d \
./src/agents/agent.d 


# Each subdirectory must supply rules for building sources it contributes
src/agents/%.o: ../src/agents/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


