# --- CONFIGURATION ---
FC = gfortran
MPIFC = mpif90
FLAGS_COMMON = -std=f2008 -Wall -O2 -J. -cpp

# specific flags:
FLAGS_SERIAL = $(FLAGS_COMMON)                      # Serial
FLAGS_OPENMP = $(FLAGS_COMMON) -fopenmp -DUSE_OMP   # OMP
FLAGS_MPI    = $(FLAGS_COMMON)                      # MPI Puro
FLAGS_HYBRID = $(FLAGS_COMMON) -fopenmp -DUSE_OMP   # Hybrid (OMP + MPI)

# source files

SOURCES_COMMON = src/geometry.f90 src/particle.f90 src/io_module.f90
SOURCES_SERIAL = $(SOURCES_COMMON) src/serial_omp/barnes_hut.f90 src/serial_omp/main.f90
# --- VERSION 1 (Original) ---
SOURCES_MPI1   = $(SOURCES_COMMON) src/mpi/barnes_hut_mpi.f90 src/mpi/main_mpi.f90
EXEC_MPI1       = barnes_hut_mpi1
EXEC_HYBRID1    = barnes_hut_hybrid1
# --- VERSION 2 (Optimizada con buffers persistentes) ---
SOURCES_MPI2   = $(SOURCES_COMMON) src/mpi/barnes_hut_mpi_v2.f90 src/mpi/main_mpi_v2.f90
EXEC_MPI2       = barnes_hut_mpi2
EXEC_HYBRID2    = barnes_hut_hybrid2

# executable names for serial/omp
EXEC_SERIAL = barnes_hut_serial
EXEC_OMP    = barnes_hut_omp

# --- EXECUTION RULES ---

help:
	@echo "Available options:"
	@echo "  make all                       -> Compile all versions"
	@echo "  make serial                    -> Compile serial version"
	@echo "  make omp                       -> Compile parallel OpenMP version"
	@echo "  make mpi1                      -> Compile MPI v1"
	@echo "  make hybrid1                   -> Compile Hybrid v1"
	@echo "  make mpi2                      -> Compile MPI v2"
	@echo "  make hybrid2                   -> Compile Hybrid v2"
	@echo "  --------------------------------------------------------"
	@echo "  make run_serial N=XX           -> Run serial version"
	@echo "  make run_omp N=XX T=4          -> Run OpenMP"
	@echo "  make run_mpi1 N=XX P=4         -> Run MPI v1"
	@echo "  make run_hybrid1 N=XX P=2 T=2  -> Run Hybrid v1"
	@echo "  make run_mpi2 N=XX P=4         -> Run MPI v2"
	@echo "  make run_hybrid2 N=XX P=2 T=2  -> Run Hybrid v2"
	@echo "  --------------------------------------------------------"
	@echo "  make input                     -> Compile and run input generator"
	@echo "  make plot N=XX                 -> Run animation"
	@echo "  make clean                     -> Clean all files"

all: serial omp mpi1 mpi2 hybrid1 hybrid2

# --- SERIAL & OMP ---
serial: $(SOURCES_SERIAL)
	@echo "Compiling SERIAL..."
	$(FC) $(FLAGS_SERIAL) -o $(EXEC_SERIAL) $(SOURCES_SERIAL)

omp: $(SOURCES_SERIAL)
	@echo "Compiling OPENMP..."
	$(FC) $(FLAGS_OPENMP) -o $(EXEC_OMP) $(SOURCES_SERIAL)

# --- MPI VERSION 1 ---
mpi1: $(SOURCES_MPI1)
	@echo "Compiling MPI Pure v1..."
	$(MPIFC) $(FLAGS_MPI) -o $(EXEC_MPI1) $(SOURCES_MPI1)

hybrid1: $(SOURCES_MPI1)
	@echo "Compiling HYBRID v1..."
	$(MPIFC) $(FLAGS_HYBRID) -o $(EXEC_HYBRID1) $(SOURCES_MPI1)

# --- MPI VERSION 2 ---
mpi2: $(SOURCES_MPI2)
	@echo "Compiling MPI Pure v2 (Optimized)..."
	$(MPIFC) $(FLAGS_MPI) -o $(EXEC_MPI2) $(SOURCES_MPI2)

hybrid2: $(SOURCES_MPI2)
	@echo "Compiling HYBRID v2 (Optimized)..."
	$(MPIFC) $(FLAGS_HYBRID) -o $(EXEC_HYBRID2) $(SOURCES_MPI2)

# --- INPUT GENERATOR ---
input_gen: generate_input_disk.f90
	$(FC) $(FLAGS_SERIAL) -o generate_input_disk generate_input_disk.f90

input: input_gen
	./generate_input_disk

# --- PARAMETRIZED EXECUTION ---
N ?= 50 # Default N
P ?= 4  # Default MPI Processes
T ?= 4  # Default OpenMP Threads (Changed default to 4 for convenience)

run_serial: serial
	@if [ ! -f input_disk_$(N).dat ]; then echo "Error: input_disk_$(N).dat missing"; exit 1; fi
	./$(EXEC_SERIAL) < input_disk_$(N).dat

# --- RUN OMP (Modified to use T) ---
run_omp: omp
	@if [ ! -f input_disk_$(N).dat ]; then echo "Error: input_disk_$(N).dat missing"; exit 1; fi
	@echo "Running OpenMP with $(T) threads..."
	export OMP_NUM_THREADS=$(T) && ./$(EXEC_OMP) < input_disk_$(N).dat 

# --- RUN MPI 1 ---
run_mpi1: mpi1
	@if [ ! -f input_disk_$(N).dat ]; then echo "Error: input_disk_$(N).dat missing"; exit 1; fi
	mpirun -np $(P) ./$(EXEC_MPI1) < input_disk_$(N).dat

run_hybrid1: hybrid1
	@if [ ! -f input_disk_$(N).dat ]; then echo "Error: input_disk_$(N).dat missing"; exit 1; fi
	export OMP_NUM_THREADS=$(T) && mpirun -np $(P) ./$(EXEC_HYBRID1) < input_disk_$(N).dat

# --- RUN MPI 2 ---
run_mpi2: mpi2
	@if [ ! -f input_disk_$(N).dat ]; then echo "Error: input_disk_$(N).dat missing"; exit 1; fi
	mpirun -np $(P) ./$(EXEC_MPI2) < input_disk_$(N).dat

run_hybrid2: hybrid2
	@if [ ! -f input_disk_$(N).dat ]; then echo "Error: input_disk_$(N).dat missing"; exit 1; fi
	export OMP_NUM_THREADS=$(T) && mpirun -np $(P) ./$(EXEC_HYBRID2) < input_disk_$(N).dat

# --- ANIMATION ---
plot:
	python3 animation.py $(N)

# --- CLEANUP ---
clean:
	rm -f *.o *.mod generate_input_disk
	rm -f $(EXEC_SERIAL) $(EXEC_OMP) 
	rm -f $(EXEC_MPI1) $(EXEC_HYBRID1) 
	rm -f $(EXEC_MPI2) $(EXEC_HYBRID2)

.PHONY: help all serial omp mpi1 mpi2 hybrid1 hybrid2 input run_serial run_omp run_mpi1 run_mpi2 run_hybrid1 run_hybrid2 plot clean