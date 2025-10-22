FC = gfortran
FFLAGS = -O2 -Wall -Wextra
OBJ = geometry.o particle.o ex1.o
EXE = ex1

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

geometry.o: geometry.f90
	$(FC) $(FFLAGS) -c geometry.f90


ex1.o: ex1.f90 geometry.o particle.o
	$(FC) $(FFLAGS) -c ex1.f90

particle.o: particle.f90 geometry.o
	$(FC) $(FFLAGS) -c particle.f90

clean:
	rm -f *.o *.mod $(EXE) output.dat

.PHONY: all clean