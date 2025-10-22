This program simulates the dynamical evolution of N particles using the Leapfrog integration method.
The initial positions, velocities, and masses of the particles are provided in the file input_ex.dat, which contains sample test data.

The program consists in three different files:

-A module called "geometry", where useful expressions of operations with vectors are defined (as well as specific types for points in 3D and vectors in 3D).

-Another module called "particle", in which the type "particle3d" is defined based on the point3d type from "geometry" storing the position, the vector3d type storing the velocity, and a real number storing the mass of the particle. 

-The code contained in the file "ex1" is the one that actually gives us the evolution of the particles with the Leapfrog method, reading the initial values from "input_ex.dat" and giving the results (with the format: time p1x p1y p1z ... pnx pny pnz) in "output.dat".

The compilation can be done with "make" using the provided Makefile to automatically build all modules and the executable.

Finally, two .gif files are provided, with animations (made with python) showing the evolution of the particles (both with and without trails).
