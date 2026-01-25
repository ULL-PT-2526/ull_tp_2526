# Makefile for Barnes–Hut N-body (serial, OpenMP, and MPI)
# Usage:
#   make            # builds serial and MPI
#   make serial     # builds serial executable
#   make omp        # builds OpenMP executable
#   make mpi        # builds MPI executable
#   make clean      # removes objects and executables

# Compilers
FC     ?= gfortran
MPIFC  ?= mpif90

# Flags
FFLAGS   ?= -O2 -Wall -Wextra
LDFLAGS  ?=

# Run configuration
# All setup files to run
SETUPS        := $(wildcard setups/*.dat)

# Defaults (can be overridden: make run OMP_THREADS=8 NP=8)
OMP_THREADS  ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
NP           ?= $(OMP_THREADS)
TIME         ?= /usr/bin/time -f %e
OUTDIR       ?= outputs
TIMEDIR      ?= times

# Sources
MOD_SRCS      := geometry.f90 particle.f90 barnes_hut_module.f90
SERIAL_MAIN   := main_barnes_hut.f90
MPI_MAIN      := main_barnes_hut_mpi.f90

# Objects (serial)
MOD_OBJS      := geometry.o particle.o barnes_hut_module.o
SERIAL_OBJS   := $(MOD_OBJS) main_barnes_hut.o
MPI_OBJS      := bh_mpi_module.o main_barnes_hut_mpi.o

# Objects (OpenMP) — compiled with -fopenmp into distinct names
OMP_MOD_OBJS  := geometry_omp.o particle_omp.o barnes_hut_module_omp.o
OMP_SERIAL_OBJS := $(OMP_MOD_OBJS) main_barnes_hut_omp.o


# Executables
SERIAL_EXE   := nbody_bh_serial
OMP_EXE      := nbody_bh_omp
MPI_EXE      := nbody_bh_mpi

.PHONY: all build serial omp mpi run run-serial run-omp run-mpi clean

all: serial

# Build everything (serial, OpenMP, MPI)
build: serial omp mpi

serial: $(SERIAL_EXE)

omp: $(OMP_EXE)

mpi: $(MPI_EXE)

# Generic rule for serial object compilation
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Serial executable
$(SERIAL_EXE): $(SERIAL_OBJS)
	$(FC) $(FFLAGS) -o $@ $(SERIAL_OBJS) $(LDFLAGS)

# Module dependency order
particle.o: geometry.o
barnes_hut_module.o: geometry.o particle.o
main_barnes_hut.o: $(MOD_OBJS)

# MPI objects: compile the MPI main 
bh_mpi_module.o: barnes_hut_mpi_module.f90
	$(MPIFC) $(FFLAGS) -c $< -o $@

main_barnes_hut_mpi.o: main_barnes_hut_mpi.f90
	$(MPIFC) $(FFLAGS) -c $< -o $@

# MPI executable
$(MPI_EXE): $(MPI_OBJS)
	$(MPIFC) $(FFLAGS) -o $@ $(MPI_OBJS) $(LDFLAGS)

# -----------------------
# OpenMP build (separate objects compiled with -fopenmp)
# -----------------------
geometry_omp.o: geometry.f90
	$(FC) $(FFLAGS) -fopenmp -c $< -o $@

particle_omp.o: particle.f90 geometry_omp.o
	$(FC) $(FFLAGS) -fopenmp -c $< -o $@

barnes_hut_module_omp.o: barnes_hut_module.f90 geometry_omp.o particle_omp.o
	$(FC) $(FFLAGS) -fopenmp -c $< -o $@

main_barnes_hut_omp.o: main_barnes_hut.f90 $(OMP_MOD_OBJS)
	$(FC) $(FFLAGS) -fopenmp -c $< -o $@

$(OMP_EXE): $(OMP_SERIAL_OBJS)
	$(FC) $(FFLAGS) -fopenmp -o $@ $(OMP_SERIAL_OBJS) $(LDFLAGS)


# Clean
clean:
	$(RM) *.o *.mod $(SERIAL_EXE) $(OMP_EXE) $(MPI_EXE) output.dat

# -----------------------
# Run all setups and store outputs/timings
# -----------------------

# Run everything (build and execute all modes over all setups)
run: build run-serial run-omp run-mpi

# Serial runs
run-serial: $(SERIAL_EXE)
	@mkdir -p $(OUTDIR) $(TIMEDIR)
	@for s in $(SETUPS); do \
	  name=$$(basename $$s .dat); \
	  name=$${name#setup_}; \
	  echo "==> [serial] $$name"; \
	  rm -f output.dat; \
	  $(TIME) -o $(TIMEDIR)/time_$${name}_serial.txt ./$(SERIAL_EXE) < $$s; \
	  mv -f output.dat $(OUTDIR)/output_$${name}_serial.dat; \
	done

# OpenMP runs
run-omp: $(OMP_EXE)
	@mkdir -p $(OUTDIR) $(TIMEDIR)
	@for s in $(SETUPS); do \
	  name=$$(basename $$s .dat); \
	  name=$${name#setup_}; \
	  echo "==> [omp $(OMP_THREADS)] $$name"; \
	  rm -f output.dat; \
	  OMP_NUM_THREADS=$(OMP_THREADS) $(TIME) -o $(TIMEDIR)/time_$${name}_omp.txt ./$(OMP_EXE) < $$s; \
	  mv -f output.dat $(OUTDIR)/output_$${name}_omp.dat; \
	done

# MPI runs
run-mpi: $(MPI_EXE)
	@mkdir -p $(OUTDIR) $(TIMEDIR)
	@for s in $(SETUPS); do \
	  name=$$(basename $$s .dat); \
	  name=$${name#setup_}; \
	  echo "==> [mpi $(NP)] $$name"; \
	  rm -f output.dat; \
	  $(TIME) -o $(TIMEDIR)/time_$${name}_mpi.txt mpirun -np $(NP) ./$(MPI_EXE) < $$s; \
	  mv -f output.dat $(OUTDIR)/output_$${name}_mpi.dat; \
	done
