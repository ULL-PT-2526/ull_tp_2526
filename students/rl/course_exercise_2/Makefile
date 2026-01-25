## Compilers
FC = gfortran
MPIFC   = mpif90

## Flags
FCFLAGS = -O2
OMPFLAGS = -fopenmp 

## Source files
#OBJECTS = geometry.o particle.o bh.o main.o
SRC_COMMON = src/geometry.f90 src/particle.f90 src/bh.f90

SRC_OMP = $(SRC_COMMON) openmp/main.f90
SRC_MPI = $(SRC_COMMON) mpi/main_mpi.f90

## Targets (default: serial)
all: nbody_bh

# In serial, I enabled openmp just to avoid issues with 'omp_get_wtime':
nbody_bh: $(SRC)
	$(FC) $(FCFLAGS) $(OMPFLAGS) -o $@ $(SRC_OMP) 

# OpenMP version
nbody_bh_omp: $(SRC)
	$(FC) $(FCFLAGS) $(OMPFLAGS) -o $@ $(SRC_OMP)

# MPI version
mpi: nbody_bh_mpi

nbody_bh_mpi: $(SRC_MPI)
	$(MPIFC) $(FCFLAGS) -o $@ $(SRC_MPI)

# Remove compiled files
clean:
	rm -f *.o *.mod nbody_bh nbody_bh_omp nbody_bh_mpi output.dat

# Run tests
test: nbody_bh
	./nbody_bh < input.dat

test_omp: nbody_bh_omp
	OMP_NUM_THREADS=4 ./nbody_bh_omp < input.dat

test_mpi: nbody_bh_mpi
	mpirun -np 4 ./nbody_bh_mpi < input.dat
