################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../libnma_affine.cpp \
../libnma_cg.cpp \
../libnma_def.cpp \
../libnma_deriv.cpp \
../libnma_diag.cpp \
../libnma_hessian.cpp \
../libnma_io.cpp \
../libnma_kinetic.cpp \
../libnma_matrix.cpp \
../libnma_misc.cpp \
../libnma_move.cpp \
../libnma_time.cpp 

OBJS += \
./libnma_affine.o \
./libnma_cg.o \
./libnma_def.o \
./libnma_deriv.o \
./libnma_diag.o \
./libnma_hessian.o \
./libnma_io.o \
./libnma_kinetic.o \
./libnma_matrix.o \
./libnma_misc.o \
./libnma_move.o \
./libnma_time.o 

CPP_DEPS += \
./libnma_affine.d \
./libnma_cg.d \
./libnma_def.d \
./libnma_deriv.d \
./libnma_diag.d \
./libnma_hessian.d \
./libnma_io.d \
./libnma_kinetic.d \
./libnma_matrix.d \
./libnma_misc.d \
./libnma_move.d \
./libnma_time.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpx -g -O3 -inline-level=2 -I../include -I../../../src -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


