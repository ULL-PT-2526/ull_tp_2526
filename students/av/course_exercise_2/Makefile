# Makefile for Course Exercise 2

#--- Compiler and Flags ---#
# Define the gfortran as Fortran Compiler (FC)
FC = gfortran

# Define general flags (debug, warnings)
FFLAGS = -O3 -Wall -fcheck=all -mmacosx-version-min=11.1

# Flags for openmp
OMPFLAGS = -fopenmp

#--- Source and Object Files ---#
# Define the list of object files (.o) explicitly
OBJS = geometry.o particle.o barnes_hut_module.o ex2.o

# Define the name of the final program
TARGET = ex2

#--- Build Rules ---#
# All depends on the target executable
all: $(TARGET)

# Rule to link the final executable
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OMPFLAGS) -o $@ $(OBJS)

# Pattern rule to compile source files into object files
%.o: %.f90
	$(FC) $(FFLAGS) $(OMPFLAGS) -c $< -o $@

#--- Module Dependencies ---#
# Dependencies between files
particle.o: geometry.o
barnes_hut_module.o: geometry.o particle.o
ex2.o: geometry.o particle.o barnes_hut_module.o

#--- Cleaning Rule ---#
# Rule to remove generated files
clean:
	rm -f $(TARGET) $(OBJS) *.mod output.dat

# Make Targets that are not actual files 
.PHONY: all clean