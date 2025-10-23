# Exercise 1: N-body simulation with leapfrog integration
This project consists of the numerical integration of a system of gravitationally interacting particles using the Leapfrog method. Its goal is to allow the simulation and visualization of the trajectories of several bodies (exact amount can be chosen by the user) that interact through gravity. Specifically, the numerical computation was written in Fortran 90 and the visualization code was written in Python.
<p align="center">
  <img src="https://github.com/user-attachments/assets/0fc81e24-00f5-466d-950b-761e5e09276d" width="350" style="margin:10px;">
  <img src="https://github.com/user-attachments/assets/14d541f8-7ea9-4648-a31e-34870674a295" width="360" style="margin:10px;">
  <img src="https://github.com/user-attachments/assets/01c5d884-1a5e-4872-af3b-6e26b5aa9e0a" width="300" style="margin:10px;">
</p>

## Project structure
```
ull_tp_2526/students/rl/course_exercise_1/
├── examples # Orbit visualizations
│   ├── simulation_broucke.gif # Broucke-Hénon orbit
│   └── simulation_circular.gif # Two bodies orbiting a common barycenter
│   └── simulation_figure8.gif # Figure 8 (Chenciner-Montgomery) orbit
├── ics # Initial conditions
│   ├── input.dat # Figure 8 (Chenciner-Montgomery) orbit
│   ├── input_broucke.txt # Broucke-Hénon orbit
│   └── input_circular.txt # Two bodies orbiting a common barycenter
│   └── input_template.txt # Template for user
├── Makefile 
├── README.md
├── geometry.f90 # Module with 3D geometry definitions and vector operations
├── particle.f90 # Module defining the particle3d type
├── ex1.f90 # Main leapfrog integration program
├── ex1_old.f90 # old version of the program, without subroutines
└── plot.py # Python script for visualizing the simulation results
```

## Code explanation

### **`geometry.f90`**

Defines the basic geometric operations that are used throughout the simulation.

Contains:
- data types vector3d and point 3d
- functions sumvp, sumpv, sumvv, sumpp, subvp, subpv, subpp, mulrv, mulvr, divvr
- operators matching these functions
- additional functions: distance, angle, norm, normalize, cross_product, dot_prod

### **`particle.f90`**

Defines a module 'particle' that contains a type 'particle3d' to represent a particle in 3D space. 
Said type has the following components:
- a point3d variable 'p' storing the particle's position
- a vector3d variable 'v' storing the particle's velocity
- a real variable 'm' storing the particle's mass  

It is used to store and manipulate particles during integration.

### **`ex1.f90`**

The main program of the simulation performs time integration through the Leapfrom method. It is set up as follows:

1. Reads the initial conditions of the system from a text file that must be passed as a command when executing the program.
2. Computes the gravitational accelerations between every pair of particles.
3. Integrates positions and velocities using the Leapfrog scheme:
4. Writes results to 'output.dat' in the format:

```
time p1x p1y p1z p2x p2y p2z ... pnx pny pnz
```
### **`plot.py`**

Python script to visualize the simulation's output in 2D. It reads the 'output.dat' file and produces an animated GIF of the system's evolution. You can install the necessary dependencies using:

pip install numpy matplotlib Pillow

## Tutorial

1. If you don't have them, install the necessary python libraries 
```
pip install numpy matplotlib pillow
```
2. Compile the program using the Makefile
```
make
```
3.  Prepare an input file (or choose one from the 'ics' folder)

You can copy the template in 'input_template.txt' and fill in your own desired values for the initial conditions of the system. Or you can also copy the format directly from here:

```
dt
dt_out
t_end
n
m1 x1 y1 z1 vx1 vy1 vz1
m2 x2 y2 z2 vx2 vy2 vz2
...
mn xn yn zn vxn vyn vzn
```

Once you have written your input file, it is advised to save it in the 'ics' folder. 

4. Run the simulation
```
./ex1 ics/input_file.txt
```
5. Visualize the results 
```
python plot.py
```
Once you have done that, the animation will be saved as 'simulation.gif'.