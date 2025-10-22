
FC = gfortran

SRCS = geometry.f90 particle.f90 ex1.f90
OBJS = geometry.o particle.o ex1.o

EXEC = ex1

all: $(EXEC)

$(EXEC): $(OBJS)
	$(FC) -o $(EXEC) $(OBJS)

geometry.o: geometry.f90
	$(FC) -c geometry.f90

particle.o: particle.f90 geometry.o
	$(FC) -c particle.f90

ex1.o: ex1.f90 geometry.o particle.o
	$(FC) -c ex1.f90

