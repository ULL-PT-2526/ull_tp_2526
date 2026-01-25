#=======================================================
# Compilers and Flags
#=======================================================
FC      = gfortran
MPIFC   = mpifort

# -cpp: Enable preprocessor (needed for #ifdef _OPENMP)
FLAGS_COMMON = -O3 -Wall -fimplicit-none -cpp
FLAGS_OPENMP = $(FLAGS_COMMON) -fopenmp
FLAGS_MPI    = $(FLAGS_COMMON)

#=======================================================
# Executables
#=======================================================
SERIAL_EXE  = nbody_serial
OMP_EXE     = nbody_openmp
MPI_EXE     = nbody_mpi

#=======================================================
# Objects
#=======================================================
MOD_OBJS    = geometry.o particle.o
BH_SERIAL_O = barnes_hut.o
BH_MPI_O    = barnes_hut_mpi.o

#=======================================================
# Run Configuration (Default: 2 Cores)
#=======================================================
N ?= 100
P ?= 2          # Procesos MPI
T ?= 2          # Hilos OpenMP

#=======================================================
# Main Rules
#=======================================================
all: serial openmp mpi

#--- SERIAL version ---
serial: $(SERIAL_EXE)

$(SERIAL_EXE): $(MOD_OBJS) $(BH_SERIAL_O) main_serial.o
	$(FC) $(FLAGS_COMMON) -o $@ $^

main_serial.o: main.f90
	$(FC) $(FLAGS_COMMON) -c main.f90 -o main_serial.o

#--- OPENMP version ---
openmp: $(OMP_EXE)

$(OMP_EXE): $(MOD_OBJS) $(BH_SERIAL_O) main_omp.o
	$(FC) $(FLAGS_OPENMP) -o $@ $^

main_omp.o: main.f90
	$(FC) $(FLAGS_OPENMP) -c main.f90 -o main_omp.o

#--- MPI version ---
mpi: $(MPI_EXE)

$(MPI_EXE): $(MOD_OBJS) $(BH_MPI_O) main_mpi.o
	$(MPIFC) $(FLAGS_MPI) -o $@ $^

main_mpi.o: main_mpi.f90
	$(MPIFC) $(FLAGS_MPI) -c main_mpi.f90 -o main_mpi.o

#=======================================================
# Module Rules (Dependencies)
#=======================================================
geometry.o: geometry.f90
	$(FC) $(FLAGS_COMMON) -c geometry.f90

particle.o: particle.f90 geometry.o
	$(FC) $(FLAGS_COMMON) -c particle.f90

barnes_hut.o: barnes_hut.f90 geometry.o particle.o
	$(FC) $(FLAGS_COMMON) -c barnes_hut.f90

barnes_hut_mpi.o: barnes_hut_mpi.f90 geometry.o particle.o
	$(MPIFC) $(FLAGS_MPI) -c barnes_hut_mpi.f90

#=======================================================
# Execution (Automated)
#=======================================================
run_serial: serial
	@echo "Running Serial Version with setup_$(N).dat..."
	@echo "'setups/setup_$(N).dat'" | ./$(SERIAL_EXE)

run_openmp: openmp
	@echo "Running OpenMP Version (T=$(T)) with setup_$(N).dat..."
	@export OMP_NUM_THREADS=$(T) && echo "'setups/setup_$(N).dat'" | ./$(OMP_EXE)

run_mpi: mpi
	@echo "Running MPI Version (P=$(P)) with setup_$(N).dat..."
	@echo "'setups/setup_$(N).dat'" | mpirun -np $(P) ./$(MPI_EXE)

#=======================================================
# Utilities
#=======================================================
clean:
	@echo "Cleaning object files and executables..."
	rm -f *.o *.mod $(SERIAL_EXE) $(OMP_EXE) $(MPI_EXE)

help:
	@echo "======================================================="
	@echo "  N-body Barnes-Hut Project "
	@echo "======================================================="
	@echo "Compilation Options:"
	@echo "  make serial       -> Compile serial version"
	@echo "  make openmp       -> Compile OpenMP version"
	@echo "  make mpi          -> Compile MPI version"
	@echo "  make all          -> Compile all versions"
	@echo ""
	@echo "Execution Options (Automated):"
	@echo "  make run_serial N=10      -> Run serial using setups/setup_10.dat"
	@echo "  make run_openmp N=100 T=2 -> Run OpenMP (100 particles, 2 threads)"
	@echo "  make run_mpi N=500 P=2    -> Run MPI (500 particles, 2 processes)"
	@echo ""
	@echo "Notes:"
	@echo "  - If N is not specified, the default is 1000."
	@echo "  - Files must be located in the 'setups/' folder."
	@echo ""
	@echo "Utilities:"
	@echo "  make clean        -> Remove executables, objects, and .mod files"
	@echo "  make help         -> Show this help"
	@echo "======================================================="

.PHONY: all serial openmp mpi run_serial run_openmp run_mpi clean help