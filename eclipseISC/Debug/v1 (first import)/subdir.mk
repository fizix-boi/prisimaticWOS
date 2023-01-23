################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../v1\ (first\ import)/ics.c 

C_DEPS += \
./v1\ (first\ import)/ics.d 

OBJS += \
./v1\ (first\ import)/ics.o 


# Each subdirectory must supply rules for building sources it contributes
v1\ (first\ import)/ics.o: ../v1\ (first\ import)/ics.c v1\ (first\ import)/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"v1 (first import)/ics.d" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-v1-20--28-first-20-import-29-

clean-v1-20--28-first-20-import-29-:
	-$(RM) ./v1\ (first\ import)/ics.d ./v1\ (first\ import)/ics.o

.PHONY: clean-v1-20--28-first-20-import-29-

