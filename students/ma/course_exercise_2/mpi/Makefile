
# Project name
NAME = ex2_mpi

# Configuration settings
FC = mpif90
# Compiler Flags
# -c: Compile only
# -fallow-argument-mismatch: Fixes the MPI_BCAST type error in gfortran 10+ since we use the example 2 that have include 'mpif.h' instead of use mpi
# -O2: Optimization level 2 makes simulation faster
FCFLAGS = -c -O2 -fallow-argument-mismatch



MOD_OBJS = geometry.o particles.o tree_mpi.o
SIM_OBJS = ex2_mpi.o # main program of the simulation
GEN_OBJS = N_body_generator.o # Inputfile generator

all: ex2_mpi generator

ex2_mpi: $(MOD_OBJS) $(SIM_OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)


generator: geometry.o N_body_generator.o
	$(FC) -o $@ $^


%.o: %.f90
	$(FC) $(FCFLAGS) $<
	# With this take all the .o and compile it first
clean:
		rm -f *.o *.mod ex2_mpi generator output.dat input.dat

# Run the executable
gen:
	./generator

run: ex2_mpi
	mpirun -np 4 ./ex2_mpi input.dat 

# We put 4 processors
