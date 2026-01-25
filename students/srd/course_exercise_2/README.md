This work is an adaptation of a Barnes-Hut code from Angel de Vicente's TAP book. The Barnes-Hut code simulates the evolution of an N-body problem.

In the classical N-body problem, the number of interactions scales with the square of the number of particles N*N, meaning that a 10x increase in particles requires 100x more computing time. The Barnes-Hut algorithm minimizes this impact by introducing an Octree to divide the space hierarchically. For distant particles, the code uses the center of mass of a tree branch as an approximation, effectively reducing the computational complexity to NlogN.

Initial conditions are generated using "particulas.f90", which distributes N particles randomly inside a sphere (currently set to radius 8) with no initial velocities. This setup physically corresponds to a "cold collapse". A softening parameter, "epsilon", was included in the force calculation to avoid divisions by near-zero distances, which prevents numerical instabilities and program crashes during close particle encounters. The output is stored in "input.dat", serving as the input for the simulation programs.

The project structure is organized as follows:

- Modules "geometry" and "particle": Define 3D vector operations and the "particle3d" derived type.
- tree_serial.f90: The sequential version of the code.
- tree_OpenMP.f90: Based on "Implementation 1" from TAP, but with OpenMP It parallelizes the main force loops. This code can be compiled without the OpenMP flag to function exactly like the serial version.
- tree_OpenMP_max.f90: An optimized version that parallelizes only the most computationally intensive loops to maximize efficiency.
- tree_MPI.f90: A distributed-memory version using MPI, based on "Implementation 1" from TAP.

Note on OpenMP: Tree construction is kept serial to avoid race conditions, where different cores might attempt to modify the same tree cell simultaneously.
The MPI implementation uses a locally redundant tree strategy, where each process constructs the full Octree after a global position synchronization (MPI_ALLGATHERV). While this increases the tree management overhead as we'll discuss in the analysis, it significantly simplifies the force calculation phase by eliminating the need for inter-process communication during tree traversal.

Execution commands:
./tree_serial < input.dat
./tree_omp < input.dat
mpirun --allow-run-as-root -np 4 ./tree_mpi < input.dat

TIME CONTROL ANALYSIS:

Performance was evaluated using three components: Tree management, Force calculation, and Integration + I/O. 

1. Small scale (500 particles):
Serial execution took 7.75s (87% in forces). OpenMP reduced this to 6.38s. At this scale, the overhead of thread assignment is noticeable. MPI was the fastest at 4.63s, despite tree management rising to 39% due to individual tree construction in each process.

2. Medium scale (5000 particles):
The theoretical O(N log N) scaling factor for a 10x increase in N is ~14. We observed:
- Serial: 138.14s (x17.8 scaling).
- OpenMP: 60.48s (x9.5 scaling). Performance is better than expected as parallelization overhead becomes negligible.
- MPI: 77.88s (x16.8 scaling).

3. Large scale (50000 particles):
- MPI: 909s (x11.67 scaling vs 5000 particles). Tree management remains the bottleneck, although force calculation is efficient.
- OpenMP: 730.2s (x12.06 scaling vs 5000 particles). At this scale, even parallelized I/O begins to show benefits, though the main gain remains in force calculation.
- Serial: 2724.90s (x19.72 scaling vs 5000 particles).

The results closely align with the N log N prediction, confirming that the algorithm becomes increasingly efficient as the system size grows. Parallelization proves essential at large scales, reducing computation time by factors of 3 to 4. For this study, OpenMP is the most effective solution due to its implementation simplicity and superior performance compared to the serial case.

Finally, an animation.gif file is provided to visualize the physical evolution of the particle set during the simulation (visualization of the cold collapse).