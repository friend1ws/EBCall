PILEUP_2_BASE_SRCS += \
subscript/Pileup2Base.cpp \
subscript/InStream.cpp

COMP_BASE_SRCS += \
subscript/CompBase.cpp \
subscript/InStream.cpp

PILEUP_2_BASE_OBJS += \
subscript/Pileup2Base.o

COMP_BASE_OBJS += \
subscript/CompBase.o

# All Target
all:PileUp2Base CompBase 

# Tool invocations
PileUp2Base: 
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -Wall -o3 -o $(PILEUP_2_BASE_OBJS) $(PILEUP_2_BASE_SRCS)
	@echo 'Finished building target: $@'
	@echo ' '

CompBase: 
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -Wall -o3 -o $(COMP_BASE_OBJS) $(COMP_BASE_SRCS)
	@echo 'Finished building target: $@'
	@echo ' '
