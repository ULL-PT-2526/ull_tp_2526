Angeles Moreno Guedes -- ma/course_exercise_1

This code has four files. The main program is ex1.f90, where the leapfrog integration is made.
The ex1.f90 uses the module geometry that includes several functions to compute vector and points operations.
Particles module, that defines a derived type to describe the principal characteristics of a particle to compute its movement
And a default input.dat file with the initial contidions.

In order to use the code, is important to first compile it with the makefile, using make in the directory where the code is:

`make`

And then if you want to run it, in the terminal:

`./ex1 <inputfile name>`

The program has a default input.dat with the initial conditions. It should have the following shape:

dt    (timestep) <br> 
dt_out     (timestep to print the result)<br> 
t     (simulation time)<br> 
n      (number of particles)<br> 
m1 x1 y1 z1 vx1 vy1 vz1 (mass, initial position, initial velocity of particle1)<br> 
.<br>
.<br>
mn xn yn zn vxn vyn vzn        (mass, initial position, initial velocity of particlen)<br> 


Finally, the output will be an output.dat with the following shape:

time    p1x     p1y     p1z     ....    pnx     pny     pnz

It will appear in the same folder when the program ends.

Besides, the folder has a out.py that you can run to make an animation from the data given by output.dat

![Image](https://github.com/user-attachments/assets/644dc391-e0ab-47bc-a2ec-1a65e1598d31)

