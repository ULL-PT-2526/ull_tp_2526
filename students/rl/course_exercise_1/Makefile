FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o ex1.o

all: ex1

ex1: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

# Remove compiled files
clean:
	rm -f *.o *.mod

# Run test
test: ex1
	./ex1 input.dat