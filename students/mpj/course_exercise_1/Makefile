# Compilator
FC = gfortran

# Options for compilation
FFLAGS = -Wall -Wextra -O2

# Files
MODULES = geometry particle
PROGRAM = ex1
OBJECTS = $(MODULES:%=%.o) $(PROGRAM).o
EXEC = ex1

# Input/output files
INPUT = input.dat
OUTPUT = output.dat

# Default rule
all: $(EXEC) run

# Compiling modules
geometry.o: geometry.f90
	$(FC) $(FFLAGS) -c $<

particle.o: particle.f90 geometry.o
	$(FC) $(FFLAGS) -c $<

# Compiling the executable
$(EXEC): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)

# Compiling the modules and main
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Exectuing the program with data
run: $(EXEC) $(INPUT)
	@echo ">>> Ejecutando $(EXEC) con $(INPUT)..."
	$(EXEC) < $(INPUT) > $(OUTPUT)
	@echo ">>> Resultado guardado en $(OUTPUT)"

# Cleaning intermediate files
clean:
	rm -f *.o *.mod $(EXEC) $(OUTPUT)
