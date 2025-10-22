This git program uses direct numerical integration to study the kinematics of a system consisting of gravitationally interacting particles using the leapfrog integration method. The principal program ex1.f90, uses the Leapfrog integration using the input.dat that has the given form:

dt

dt_out

t_end

n

m1 p1x p1y p1z v1x v1y v1z

m2 p2x p2y p2z v2x v2y v2z

...

mn pnx pny pnz vnx vny vnz

where:

- **`dt`** = Integration time step  
- **`dt_out`** = Time interval for data output  
- **`t_end`** = Total simulation time  
- **`n`** = Number of particles in the system  
- **`m`** = Particle mass  
- **`p(x, y, z)`** = Initial position vector of each particle  
- **`v(x, y, z)`** = Initial velocity vector of each particle 

The program ex1 also uses two modules: the geometry and particle module.

The geometry module "geometry.f90" contains the definition of types vector3d and point3d, both with (64-bit) real components x, y, and z the functions sumvp, sumpv, subvp, subpv, mulrv, mulvr, divvr for adding and subtracting points and vectors, and for multiplying and dividing vectors with reals. Also the functions of the distance between two points, the angle between two vectors, the nomalization of a vector and the cross product of two vectors.

The particle module "particle.f90" also uses geometry.f90 and define a type particle3d. This type have the components: point3d variable p (storing the particle’s position), vector3d   variable v (storing the particle’s velocity), and a real variable m (storing the particle’s mass).

To use the code is only necesary puting command in the terminal to compile:

make

Then to run the code and get the output:

./ex1 input.dat

I have an input data of example in the fold already, but to use other input it is only necesary to replace the following of ./ex1 for the name and location of the new file.

As said, the result will be a an output.dat file with the time positions and velocities in all the (x,y,z) directions. Attached to the fold there is also an animation showing the particles trajectories.
