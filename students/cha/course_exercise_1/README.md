# N-BODY SIMULATION USING LEAPFROG INTEGRATION IN *FORTRAN*
<span style="font-size:1.5em;">Álvaro Cascales Hernández</span>

This project implements a gravitational simulation of the well-known N-body problem. The numerical core is written in modern *Fortran* and uses a Leapfrog integrator to compute the temporal evolution of the particles. In addition, a *Python* script is included to generate a 3D animated visualization of the resulting trajectories.

## Main Features:

- **Fortran core:** high efficiency for numerical calculations in this type of simulation.  
- **Modular design:** the code is organized into modules (`geometry.f90`, `particle.f90`) to improve structure and reusability.  
- **Leapfrog integration:** a symplectic, second-order numerical method, ideal for orbital mechanics systems due to its excellent energy conservation properties over time.  
- **Softening parameter:** includes a constant $\epsilon$ to prevent numerical divergences caused by singularities when two particles get too close. This parameter is set to 0 in our case, as the input file does not cause any issues after visualization. It is recommended not to enable this parameter unless strictly necessary.  
- **3D visualization:** a *Python* script using `matplotlib` and `FuncAnimation` generates a video (`.mp4`) that shows the system’s dynamics in three dimensions, including particle trails for better visualization.  

## Project Structure:

The project is organized as follows:
- `geometry.f90`: module defining datatypes for `vector3d` and `point3d`, and overloading operators (+, -, *, /) to allow intuitive mathematical syntax.  
- `particle.f90`: module defining the `particle3d` type, which encapsulates a particle’s mass, position, and velocity in three dimensions.  
- `ex1.f90`: the first functional version of the main program. It generates an `output.dat` file.  
- `ex1v2.f90`: an improved version of the main program, including several changes inspired by classmates’ feedback and personal refinements.  
- `Makefile`: simplifies compilation and cleanup of generated files.  
- `animation.py`: *Python* script that reads the generated `output.dat` file and produces the 3D animation.  
- `input.dat`: example file containing the initial conditions for the simulation.  
- `output.dat`: output file containing the position of each particle at every time step.  

## How to Use the Simulator

Follow these steps to compile and run the project from the terminal:

1. **Prerequisites**  
Ensure you have a *Fortran* compiler (*gfortran*) and *Python 3* installed. You will also need the following *Python* libraries: `numpy`, `matplotlib`, and `ffmpeg`.

2. **Compilation**  
The included `Makefile` simplifies the process. To compile the `ex1v2` program (the recommended version), simply run:

```bash
make
```

This will create the executable file `ex1v2`.  
If you wish to clean up the generated object (`.o`) and module (`.mod`) files, run:

```bash
make clean
```

3. **Running the Simulation**  
The main program expects the name of the input file containing the initial conditions as an argument. To run the simulation using the example file `input.dat`:

```bash
./ex1v2 input.dat
```

The program will start computing and save the results to the `output.dat` file.  
The format of the input file is as follows:

```text
<dt>
<dt_out>
<t_end>
<n>
<m1> <x1> <y1> <z1> <vx1> <vy1> <vz1>
<m2> <x2> <y2> <z2> <vx2> <vy2> <vz2>
...
```

4. **Viewing Results**  
Once the simulation finishes and `output.dat` has been generated, you can view the results by running the *Python* script:

```bash
python3 animation.py
```

The script will ask if you want to save the animation as a video.  
If you answer “yes”, an `animation.mp4` file will be generated in the same directory.  
Note that saving the animation as a video (instead of only displaying it) may take longer.  

## Code Details
### Fortran Modules (`geometry.f90` and `particle.f90`)
The code is structured modularly to promote clarity and reusability:

- The `geometry` module handles all vector operations, overloading operators to allow quasi-mathematical notation (e.g., `vector_a + vector_b`). It also includes functions such as `distance2` for efficient squared-distance calculation.  
- The `particle` module relies on `geometry` to define a particle with its fundamental physical properties.  

### Main Program (`ex1v2.f90`)
The `ex1v2` version refactors the initial version, improving code structure and readability.  
The main enhancement is the `compute_accelerations` subroutine, which isolates the force calculation logic, making the main loop more readable.  
This subroutine computes accelerations for all particles based on Newton’s law of universal gravitation.  
Other minor improvements are also included compared to the first version.  

### Leapfrog Integrator
The main loop implements the *Leapfrog* “kick-drift-kick” integration scheme:

1. *Kick:* update velocities by half a time step  
   $\vec{v}_{i+1/2} = \vec{v}_i + \vec{a}_i \frac{\Delta t}{2}$  
2. *Drift:* update positions using the new velocities  
   $\vec{p}_{i+1} = \vec{p}_i + \vec{v}_{i+1/2} \Delta t$  
3. *Acceleration calculation:* compute the new accelerations $\vec{a}_{i+1}$ at the new positions $\vec{p}_{i+1}$  
4. *Kick:* complete the velocity update  
   $\vec{v}_{i+1} = \vec{v}_{i+1/2} + \vec{a}_{i+1} \frac{\Delta t}{2}$

