FC = gfortran
FFLAGS = -O2 -Wall -std=f2008
OBJ = geometry.o particle.o ex1.o
EXEC = ex1

all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

geometry.o: geometry.f90
	$(FC) $(FFLAGS) -c geometry.f90

particle.o: particle.f90 geometry.o
	$(FC) $(FFLAGS) -c particle.f90

ex1.o: ex1.f90 particle.o geometry.o
	$(FC) $(FFLAGS) -c ex1.f90

clean:
	rm -f *.o *.mod $(EXEC)

mrproper: clean
	rm -f output.dat