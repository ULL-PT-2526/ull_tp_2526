# Fortran lang has a tutorial to create a makefile
# Project name
NAME = ex2

# Configuration settings
FC = gfortran
FCFLAGS = -c -O3 -cpp

# To use OpenMP we can write USE_OMP=1. Or in terminal: make USE_OMP=1
USE_OMP ?= 0
OMPFLAG = -fopenmp

ifeq ($(USE_OMP),1)
    FCFLAGS += $(OMPFLAG)
    LDFLAGS = $(OMPFLAG)
else
    LDFLAGS =
endif

MOD_OBJS = geometry.o particles.o tree.o
SIM_OBJS = ex2.o # main program of the simulation
GEN_OBJS = N_body_generator.o # Inputfile generator

all: ex2 generator

ex2: $(MOD_OBJS) $(SIM_OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)


generator: geometry.o N_body_generator.o
	$(FC) -o $@ $^


%.o: %.f90
	$(FC) $(FCFLAGS) $<
	# With this take all the .o and compile it first
clean:
		rm -f *.o *.mod ex2 generator output.dat input.dat

# Run the executable
gen:
	./generator

run: ex2
	./ex2 input.dat

