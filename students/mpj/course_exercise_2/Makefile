# ==========================================
# Compilers and Flags
# ==========================================
FC_SERIE = gfortran
FC_OMP   = gfortran
FC_MPI   = mpif90
FFLAGS   = -O3 -Wall -fimplicit-none
OMP_FLAG = -fopenmp

# Executables
GEN_EXE   = input_rotation
EXE_SERIE = tree_serie
EXE_OMP   = tree_omp
EXE_MPI   = tree_mpi

# ==========================================
# Main Rules
# ==========================================
all: $(GEN_EXE) $(EXE_SERIE) $(EXE_OMP) $(EXE_MPI)

# 1. Data generator
$(GEN_EXE): input_rotation.f90
	$(FC_SERIE) $(FFLAGS) -o $@ $<

# 2. Serial version
$(EXE_SERIE): geometry_s.o particle_s.o tree_serial.o
	$(FC_SERIE) $(FFLAGS) -o $@ $^

# 3. OpenMP version
$(EXE_OMP): geometry_omp.o particle_omp.o tree_openmp.o
	$(FC_OMP) $(FFLAGS) $(OMP_FLAG) -o $@ $^

# 4. MPI version
$(EXE_MPI): geometry_mpi.o particle_mpi.o tree_mpi.o
	$(FC_MPI) $(FFLAGS) -o $@ $^

# ==========================================
# Object Compilation
# ==========================================

# NOTE: When compiling the first object of each set, the .mod file is generated.
# We add the manual dependency to ensure the order.

# Serial objects
geometry_s.o: geometry.f90
	$(FC_SERIE) $(FFLAGS) -c $< -o $@
particle_s.o: particle.f90 geometry_s.o
	$(FC_SERIE) $(FFLAGS) -c $< -o $@
tree_serial.o: tree_serial.f90 particle_s.o
	$(FC_SERIE) $(FFLAGS) -c $< -o $@

# OpenMP objects
geometry_omp.o: geometry.f90
	$(FC_OMP) $(FFLAGS) $(OMP_FLAG) -c $< -o $@
particle_omp.o: particle.f90 geometry_omp.o
	$(FC_OMP) $(FFLAGS) $(OMP_FLAG) -c $< -o $@
tree_openmp.o: tree_openmp.f90 particle_omp.o
	$(FC_OMP) $(FFLAGS) $(OMP_FLAG) -c $< -o $@

# MPI objects
geometry_mpi.o: geometry.f90
	$(FC_MPI) $(FFLAGS) -c $< -o $@
particle_mpi.o: particle.f90 geometry_mpi.o
	$(FC_MPI) $(FFLAGS) -c $< -o $@
tree_mpi.o: tree_MPI.f90 particle_mpi.o
	$(FC_MPI) $(FFLAGS) -c $< -o $@

# ==========================================
# Execution and Benchmarking
# ==========================================
run_all: all
	@echo "--- Generating input_rotation.dat ---"
	@./$(GEN_EXE)
	
	@echo "\n--- Running SERIAL VERSION ---"
	@./$(EXE_SERIE) < input_rotation.dat
	
	@echo "\n--- Running OPENMP ---"
	@export OMP_NUM_THREADS=8 && ./$(EXE_OMP) < input_rotation.dat
	
	@echo "\n--- Running MPI ---"
	@mpirun -np 4 ./$(EXE_MPI) < input_rotation.dat

clean:
	rm -f *.o *.mod $(GEN_EXE) $(EXE_SERIE) $(EXE_OMP) $(EXE_MPI) input_rotation.dat output*.dat