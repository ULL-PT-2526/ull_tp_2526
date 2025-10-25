# Compiler and flags
FC = gfortran
FFLAGS = -Wall -O2

# Source and object files
SRC = geometry.f90 particle.f90 ex1.f90
OBJ = geometry.o particle.o ex1.o

# Executable name
EXEC = ex1

# Default rule
all: $(EXEC)

# How to build the executable
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

# Module dependencies
particle.o: geometry.o
ex1.o: geometry.o particle.o

# Generic compilation rule
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Rule to run the program
test: $(EXEC)
	./$(EXEC)

# Clean up
clean:
	rm -f *.o *.mod $(EXEC)
