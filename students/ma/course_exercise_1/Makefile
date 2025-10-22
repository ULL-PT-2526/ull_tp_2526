# Fortran lang has a tutorial to create a makefile
# Project name
NAME = excercise1

# Configuration settings
FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particles.o ex1.o

all: ex1

ex1: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS) 
	#With this we compile and link the objects

%.o: %.f90
	$(FC) $(FCFLAGS) $<
	# With this take all the .o and compile it first
clean:
		rm -f *.o *.mod ex1

# Run the executable
test: ex1
	./ex1

