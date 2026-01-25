FC = gfortran
MPIFC = mpif90
FFLAGS = -O3
OMPFLAGS = -fopenmp

MOD_OBJS = geometry.o particle.o

EXEC_SERIAL = tree_serial
EXEC_OMP = tree_omp
EXEC_OMP_MAX = tree_omp_max
EXEC_MPI = tree_mpi

all: $(EXEC_SERIAL) $(EXEC_OMP) $(EXEC_OMP_MAX) $(EXEC_MPI)


$(EXEC_SERIAL): $(MOD_OBJS) tree_serial.o
	$(FC) $(FFLAGS) -o $(EXEC_SERIAL) $(MOD_OBJS) tree_serial.o

$(EXEC_OMP): $(MOD_OBJS) tree_OpenMP.o
	$(FC) $(FFLAGS) $(OMPFLAGS) -o $(EXEC_OMP) $(MOD_OBJS) tree_OpenMP.o

$(EXEC_OMP_MAX): $(MOD_OBJS) tree_OpenMP_max.o
	$(FC) $(FFLAGS) $(OMPFLAGS) -o $(EXEC_OMP_MAX) $(MOD_OBJS) tree_OpenMP_max.o

$(EXEC_MPI): $(MOD_OBJS) tree_MPI.o
	$(MPIFC) $(FFLAGS) -o $(EXEC_MPI) $(MOD_OBJS) tree_MPI.o

geometry.o: geometry.f90
	$(FC) $(FFLAGS) -c geometry.f90

particle.o: particle.f90 geometry.o
	$(FC) $(FFLAGS) -c particle.f90

tree_serial.o: tree_serial.f90 geometry.o particle.o
	$(FC) $(FFLAGS) -c tree_serial.f90

tree_OpenMP.o: tree_OpenMP.f90 geometry.o particle.o
	$(FC) $(FFLAGS) $(OMPFLAGS) -c tree_OpenMP.f90

tree_OpenMP_max.o: tree_OpenMP_max.f90 geometry.o particle.o
	$(FC) $(FFLAGS) $(OMPFLAGS) -c tree_OpenMP_max.f90

tree_MPI.o: tree_MPI.f90 geometry.o particle.o
	$(MPIFC) $(FFLAGS) -c tree_MPI.f90

clean:
	rm -f *.o *.mod $(EXEC_SERIAL) $(EXEC_OMP) $(EXEC_OMP_MAX) $(EXEC_MPI)
