## =========================================
# Makefile - Barnes-Hut N-body Simulation
# =========================================

# Compilers
FC = gfortran
MPIFC = mpif90

# Compiler flags
FFLAGS = -O3 -Wall

# Source files
SOURCES = geometry.f90 particle.f90 barnes_hut.f90 main.f90
MPI_SOURCES = geometry.f90 particle.f90 barnes_hut.f90 main_mpi.f90

# Executables
SERIAL = nbody_serial
PARALLEL = nbody_parallel  
MPI = nbody_mpi

# Default: compile all versions
all: serial parallel mpi

# Serial version (no OpenMP)
serial:
	$(FC) $(FFLAGS) -o $(SERIAL) $(SOURCES)
	@echo "✓ Compiled: $(SERIAL)"

# Parallel version (with OpenMP)
parallel:
	$(FC) $(FFLAGS) -fopenmp -o $(PARALLEL) $(SOURCES)
	@echo "✓ Compiled: $(PARALLEL)"

# MPI version
mpi:
	$(MPIFC) $(FFLAGS) -o $(MPI) $(MPI_SOURCES)
	@echo "✓ Compiled: $(MPI)"

# Clean
clean:
	rm -f *.o *.mod $(SERIAL) $(PARALLEL) $(MPI) output*.dat
	@echo "✓ Cleaned"

# Help
help:
	@echo "Usage:"
	@echo "  make           - Compile all versions"
	@echo "  make serial    - Compile serial version"
	@echo "  make parallel  - Compile OpenMP version"
	@echo "  make mpi       - Compile MPI version"
	@echo "  make clean     - Remove compiled files"
	@echo ""
	@echo "To run and measure time:"
	@echo "  time ./nbody_serial < input.dat"
	@echo "  export OMP_NUM_THREADS=4 && time ./nbody_parallel < input.dat"
	@echo "  time mpirun --allow-run-as-root -np 4 ./nbody_mpi < input.dat"

.PHONY: all serial parallel mpi clean help