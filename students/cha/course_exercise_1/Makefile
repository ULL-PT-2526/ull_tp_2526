# Define the compiler and flags
FC = gfortran
FFLAGS = -Wall -O2

# Default target: builds the final executable
all: ex1v2

# Rule to compile the geometry module
geometry.o: geometry.f90
	$(FC) $(FFLAGS) -c geometry.f90

# Rule to compile the particle module
# Depends on geometry module too
particle.o: particle.f90 geometry.o
	$(FC) $(FFLAGS) -c particle.f90

# Rule to compile the main program object file
ex1v2.o: ex1v2.f90 geometry.o particle.o
	$(FC) $(FFLAGS) -c ex1v2.f90

# Rule to link all object files into the final executable
ex1v2: ex1v2.o geometry.o particle.o
	$(FC) $(FFLAGS) -o ex1v2 ex1v2.o geometry.o particle.o

# Rule to clean up the directory
clean:
	rm -f *.o *.mod ex1v2




