# Barnes-Hut N-body Simulation

An implementation of the Barnes–Hut algorithm for efficient gravitational N-body simulations with serial, parallel (OpenMP), and distributed (MPI) implementations.

## Project Structure

```
├── Makefile              # Build system
├── README.md             # This file
├── geometry.f90          # 3D vector/point operations
├── particle.f90          # Particle type definition
├── barnes_hut.f90        # Barnes-Hut tree algorithm
├── main.f90              # Serial/OpenMP main program
├── main_mpi.f90          # MPI version
├── comparison_times.py  # Input file generator and comparation times
```

## Algorithm Overview

The Barnes-Hut algorithm reduces N-body force calculation complexity from O(N²) to O(N log N) using an octree:

1. **Tree Construction**: Divides recursively 3D space into octants, creating a tree where each leaf contains ≤1 particle
2. **Mass Calculation**: Calculate total mass and center of mass for each tree node
3. **Force Calculation**: For each particle, traverse the tree:
   - If a node is far enough (l/D < θ), treat it as a single mass
   - Otherwise, recurse into its children
4. **Time Integration**: Update positions and velocities using leapfrog integration

- **θ (theta) = 1.0**: Opening angle criterion (smaller = more accurate but slower)


## Compilation

### Build Commands

```bash
# Compile all versions
make

# Compile specific versions
make serial      # Serial version
make parallel    # OpenMP version
make mpi         # MPI version

# Clean compiled files
make clean

# Show help
make help
```

## Execution

### 1. Serial Version

```bash
./nbody_serial < input.dat
```

### 2. OpenMP Parallel Version

```bash
# Set number of threads
export OMP_NUM_THREADS=4

# Run
./nbody_parallel < input.dat
```

To test with different thread counts:
```bash
export OMP_NUM_THREADS=1 && ./nbody_parallel < input.dat
export OMP_NUM_THREADS=2 && ./nbody_parallel < input.dat
export OMP_NUM_THREADS=4 && ./nbody_parallel < input.dat
export OMP_NUM_THREADS=8 && ./nbody_parallel < input.dat
```

### 3. MPI Distributed Version

```bash
# Run with 4 processes
mpirun --allow-run-as-root -np 4 ./nbody_mpi < input.dat
```

To test with different process counts:
```bash
mpirun --allow-run-as-root -np 1 ./nbody_mpi < input.dat
mpirun --allow-run-as-root -np 2 ./nbody_mpi < input.dat
mpirun --allow-run-as-root -np 4 ./nbody_mpi < input.dat
mpirun --allow-run-as-root -np 8 ./nbody_mpi < input.dat
```
Note: `--allow-run-as-root` is only needed when running MPI as root (e.g. in WSL or containers).

Standard execution:

```bash
mpirun -np 4 ./nbody_mpi < input.dat
```

## Creating Input Files

### Input File Format


```bash
dt          # Time step 
dt_out      # Output interval 
t_end       # End time 
n           # Number of particles 
# For each particle:
mass x y z vx vy vz
```

### Using the Input Generator and benchmark of times

Test all versions with multiple problem sizes:

```bash
python3 comparison_times.py
```

This will:
- Generate input files with 100, 1000, 5000, 10000  particles
- Run serial, OpenMP (2, 4, 8 threads), and MPI (2, 4, 8 processes) versions
- Calculate speedups and efficiencies
- Display results table

## Output Files

- `output.dat`: Serial/OpenMP results
- `output_mpi.dat`: MPI results

Format: Each line contains `time x1 y1 z1 x2 y2 z2 ... xn yn zn`

## Parallelization Strategies

### OpenMP (Shared Memory)
- **Parallelized section**: the force calculation step 

Each thread is responsible for computing the forces acting on a group of particles.

Synchronization is not required, since each particle’s acceleration is computed independently.

```fortran
!$omp parallel do default(shared) private(i) schedule(dynamic)
do i = 1, n
    call calculate_forces_aux(i, head, particles, accelerations)
end do
!$omp end parallel do
```

### MPI (Distributed Memory)
The set of particles is divided among different processes. Every process builds its own Barnes–Hut tree, which is needed to compute forces  
- **Communication between processes**:
  - `MPI_Bcast` is used to share particle data  
  - `MPI_Allreduce` is used to combine acceleration results from all processes  

In both cases, it is necessary to take into account the overhead, an additional cost introduced by thread or process management, synchronization, and communication when a program is executed in parallel.

## Performance Results

Results obtained on a system testing different particle counts:

### Speedup Comparison

| Particle Count | Serial (s) | OpenMP 2T | OpenMP 4T | OpenMP 8T | MPI 2P | MPI 4P | MPI 8P |
|----------------|------------|-----------|-----------|-----------|--------|--------|--------|
| **N = 100**    | 0.09       | 0.77×     | 0.74×     | 1.14×     | 0.18×  | 0.16×  | 0.13×  |
| **N = 1000**   | 2.56       | 0.98×     | 1.21×     | 1.66×     | 1.04×  | 1.29×  | 1.69×  |
| **N = 5000**   | 24.62      | 0.98×     | 1.25×     | 1.82×     | 1.45×  | 2.09×  | **2.80×** |
| **N = 10000**  | 60.97      | 1.04×     | 1.31×     | 1.97×     | 1.52×  | 2.24×  | **3.03×** |

*Speedup = Serial Time / Parallel Time*

### Observations

1. **Problem Size Dependency**
   - **Small problems (N = 100)**: Parallelization overhead dominates, serial is fastest
   - **Medium problems (N = 1000)**: OpenMP and MPI become competitive (1.2-1.7× speedup)
   - **Large problems (N ≥ 5000)**: Clear benefit from parallelization
   - **Very large problems (N = 10000)**: **MPI 8P achieves 3.03× speedup** (best result)
   
2. **OpenMP Performance**
   - Consistent scaling improvement with thread count
   - N=1000: 1.66× speedup with 8 threads
   - N=10000: 1.97× speedup with 8 threads
   - Low overhead makes it effective even for medium-sized problems
   
3. **MPI Performance**
   - Large underperformance for small problems (N=100: 0.13× with 8P)
   - Becomes better at N=1000 (1.69× with 8P)
   - **Outperforms OpenMP at N≥5000**:
     - N=5000: 2.80× vs 1.82× (OpenMP 8T)
     - N=10000: 3.03× vs 1.97× (OpenMP 8T)
   - Scales better than OpenMP for large problems
   
4. **Scalability Analysis**
   - **OpenMP**: Sub-linear scaling, efficiency decreases with more threads
     - N=10000, 8T: 1.97× (24.6% efficiency)
   - **MPI**: Better strong scaling for large problems
     - N=10000, 8P: 3.03× (37.9% efficiency)
   - **Overhead**: Visible at small N but amortized at large N
   
### Recommendations

- **For N < 1000**: Use serial version (minimal overhead)
- **For N = 1000-5000**: Use OpenMP with 4-8 threads (good speedup, simple)
- **For N > 5000**: Use MPI with 8 processes (best performance)

## Technical Details

### Barnes-Hut Tree Structure
- **Cell types**: 
  - 0 = empty
  - 1 = leaf (single particle)
  - 2 = internal node (8 children)
- **Tree rebuild**: Complete rebuild every timestep after position updates
- **Memory management**: Cleanup with `borrar_tree` to prevent leaks

### Time Integration
Uses leapfrog method:
```
v(t+dt/2) = v(t) + a(t) × dt/2
r(t+dt) = r(t) + v(t+dt/2) × dt
a(t+dt) = calculate_forces(r(t+dt))
v(t+dt) = v(t+dt/2) + a(t+dt) × dt/2
```

## Author
Cristina Bolaños Quevedo  

Course Exercise 2 – Programming Techniques  
Barnes–Hut N-body Simulation  
Universidad de La Laguna, 2025–2026