# Barnes-Hut Algorithm - N-Body Simulation - Course Exercise 2

## Description

This project implements an advanced N-body gravitational simulation in Fortran using the Barnes-Hut algorithm. Unlike the direct O(N^2) approach, this implementation uses an Octree data structure to group distant particles into centers of mass, reducing the computational complexity to O(NlogN). This allows for the simulation of significantly larger systems.

The simulation supports three execution modes:

1. **Serial**: Standard sequential execution.
2. **OpenMP**: Parallel execution on shared-memory systems (multi-core CPUs).
3. **MPI**: Distributed execution for clusters or multi-process environments.

The physical integration is performed using the Velocity Verlet (Leapfrog) scheme to ensure numerical stability and energy conservation.

## Files

- `ex2.f90`: Main program for Serial and OpenMP execution. Handles command-line arguments for thread count (`-t`) and I/O.

- `ex2_mpi.f90`: Main program for the distributed MPI version. Implements domain decomposition and process synchronization via `MPI_ALLREDUCE`.

- `barnes_hut_module.f90`: The core logic of the simulation. Contains the Octree structure, recursive tree building, and the gravitational force calculation using the multipole expansion approximation.

- `geometry.f90`: Defines `vector3d` and `point3d` types and their associated mathematical operators.

- `particle.f90`: Defines the `particle3d` data structure.

- `Makefile`: Used for compiling the Serial and OpenMP versions (`ex2`).

- `Makefile_mpi`: Used for compiling the MPI version (`ex2_mpi`) using the `mpif90` wrapper.

- `ex2_animation.py`: Python script to visualize results as a 3D animation.

## Performance Comparison

The following table summarizes the execution times (in seconds) for different particle counts across the three configurations.

| N Particles   | Serial (1 thread) | OpenMP (10 threads) | MPI (10 procs) | MPI (Single ALLREDUCE) |
|---------------|-------------------|---------------------|----------------|------------------------|
| Default (3)   | 0.037806s         | 0.274575s           | 0.107451s      | 0.007487s              |
| 10            | 0.117000s         | 0.371804s           | 0.151225s      | 0.037927s              |
| 100           | 1.570561s         | 1.438253s           | 1.409588s      | 1.799579s              |
| 1000          | 48.883531s        | 18.13809s           | 39.94954s      | 93.08835s              |
| 10000         | 3220.5481s        | 604.1928s           | 1399.597s      |                        |

> **Note** Benchmarks were performed on a MacBook Pro M3 Pro with 12 cores. MPI and OpenMP times depend on the number of processes/threads used. MPI shows higher overhead in local machines due to redundant tree construction and communication latency.

### Perfomance analysis

The results demonstrate the trade-off between computational load and parallelization overhead. For small systems (N≤100), the Serial version is most efficient, as the cost of spawning threads or processes outweighs the calculation benefits. As the system size grows, parallelization becomes critical. OpenMP proves to be the optimal solution for this single-node setup, achieving a ~5.3x speedup at N=10,000. While MPI significantly outperforms the serial version (achieving a ~2.3x speedup), it shows higher latency than OpenMP in this local environment. This is due to the overhead of memory copying during `MPI_ALLREDUCE` operations and the redundant construction of the full Octree by every independent process.

### Autogenerate particle setup

To generate a specific setup, use:

```Bash
python generate_input.py <number_of_particles>
```

## Prerequisites

- **Fortran Compiler**: `gfortran` (for Serial/OpenMP) and `mpif90` (for MPI).

- **MPI Library**: OpenMPI or MPICH installed in your environment.

- **Python 3.x**: With `numpy` and `matplotlib` for visualization.

- **FFmpeg**: Required to save the `.mp4` animations.

## Compilation

The project uses two separate Makefiles to manage the different compilation requirements:

- **To compile Serial/OpenMP version**:

```Bash
make clean
make
```

- **To compile MPI version**:

```Bash
make -f Makefile_mpi clean
make -f Makefile_mpi
```

## Running the Simulation

### Serial / OpenMP

You can specify the number of threads using the `-t` flag:

```Bash
./ex2 -i input.dat -o output.dat -t 10
```

### MPI

Use `mpirun` to launch the distributed version. The number of processes is defined by the `-np` flag:

```Bash
mpirun -np 10 ./ex2_mpi -i input.dat -o output.dat
```

## How the Code Works (General Logic)

1. **Initialization**: The Master process (Rank 0 in MPI) reads the `input.dat` file containing the initial masses, positions, and velocities.

2. **Tree Construction**: In every time step, the code clears the previous Octree and builds a new one based on the current particle positions. In MPI, every process builds the full tree to avoid excessive communication during force calculation.

3. **Force Calculation**:

    - **Serial/OpenMP**: Particles are iterated, and the tree is traversed to calculate the force. OpenMP divides these iterations among threads.

    - **MPI**: Each process is assigned a subset of particles. After calculating local forces, `MPI_ALLREDUCE` is used to synchronize the global acceleration array across all processes.

4. **Integration**: Positions and velocities are updated using the Velocity Verlet algorithm.

5. **Output**: Every `dt_out` seconds, the updated positions are written to the output file (handled exclusively by Rank 0 in the MPI version).

## Visualization

Generate the 3D animation using:

```Bash
python ex2_animation.py -i output.dat -o output_simulation
```
