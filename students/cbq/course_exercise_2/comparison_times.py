import subprocess
import time
import random
import os


# Generation of input files

def generate_input(n_particles, filename):
    print(f"Generating {filename} with {n_particles} particles...")
    with open(filename, 'w') as f:
        f.write("0.01\n")   # dt
        f.write("0.5\n")    # dt_out
        f.write("5.0\n")    # t_end
        f.write(f"{n_particles}\n")
        
        for i in range(n_particles):
            mass = random.uniform(0.5, 2.0)
            x = random.uniform(-10, 10)
            y = random.uniform(-10, 10)
            z = random.uniform(-10, 10)
            vx = random.uniform(-1, 1)
            vy = random.uniform(-1, 1)
            vz = random.uniform(-1, 1)
            f.write(f"{mass} {x} {y} {z} {vx} {vy} {vz}\n")


# Running simulation and measure time

def run_simulation(command, input_file):
    with open(input_file, 'r') as f:
        start = time.perf_counter()
        subprocess.run(command, stdin=f, stdout=subprocess.DEVNULL, 
                      stderr=subprocess.DEVNULL, shell=True)
        return time.perf_counter() - start


# Main script

print("Barnes-Hut Comparison Tool")
print()

particle_counts = [100, 1000, 5000, 10000]


for n in particle_counts:
    generate_input(n, f"input_{n}.dat")

all_results = []

for n in particle_counts:
    print(f"\nTesting N = {n} particles:")
    input_file = f"input_{n}.dat"
    
    # Serial
    print(f"  Serial...     ", end="", flush=True)
    t_serial = run_simulation("./nbody_serial", input_file)
    print(f"{t_serial:.4f}s")
    
    # OpenMP 2 threads
    os.environ["OMP_NUM_THREADS"] = "2"
    print(f"  OpenMP (2T)...", end="", flush=True)
    t_omp2 = run_simulation("./nbody_parallel", input_file)
    print(f"{t_omp2:.4f}s")
    
    # OpenMP 4 threads
    os.environ["OMP_NUM_THREADS"] = "4"
    print(f"  OpenMP (4T)...", end="", flush=True)
    t_omp4 = run_simulation("./nbody_parallel", input_file)
    print(f"{t_omp4:.4f}s")

    # OpenMP 8 threads
    os.environ["OMP_NUM_THREADS"] = "8"
    print(f"  OpenMP (8T)...", end="", flush=True)
    t_omp8 = run_simulation("./nbody_parallel", input_file)
    print(f"{t_omp8:.4f}s")
    
    # MPI 2 processes
    print(f"  MPI (2P)...   ", end="", flush=True)
    t_mpi2 = run_simulation("mpirun --allow-run-as-root -np 2 ./nbody_mpi", input_file)
    print(f"{t_mpi2:.4f}s")
    
    # MPI 4 processes
    print(f"  MPI (4P)...   ", end="", flush=True)
    t_mpi4 = run_simulation("mpirun --allow-run-as-root -np 4 ./nbody_mpi", input_file)
    print(f"{t_mpi4:.4f}s")

    # MPI 8 processes
    print(f"  MPI (8P)...   ", end="", flush=True)
    t_mpi8 = run_simulation("mpirun --allow-run-as-root -np 8 ./nbody_mpi", input_file)
    print(f"{t_mpi8:.4f}s")

    all_results.append({'n': n,'serial': t_serial,'omp2': t_omp2,'omp4': t_omp4,
                        'omp8': t_omp8,'mpi2': t_mpi2,'mpi4': t_mpi4,'mpi8': t_mpi8})

# Comparison table
print()
print("=" * 60)
print("RESULTS SUMMARY")
print("=" * 60)
print()

for result in all_results:
    n = result['n']
    t_s = result['serial']
    
    print(f"N = {n} particles:")
    print(f"  {'Version':<15} {'Time':>8}  {'Speedup':>8}")
    print(f"  {'-'*15} {'-'*8}  {'-'*8}")
    print(f"  {'Serial':<15} {t_s:>7.4f}s  {1.0:>7.2f}x")
    print(f"  {'OpenMP (2T)':<15} {result['omp2']:>7.4f}s  {t_s/result['omp2']:>7.2f}x")
    print(f"  {'OpenMP (4T)':<15} {result['omp4']:>7.4f}s  {t_s/result['omp4']:>7.2f}x")
    print(f"  {'OpenMP (8T)':<15} {result['omp8']:>7.4f}s  {t_s/result['omp8']:>7.2f}x")
    print(f"  {'MPI (2P)':<15} {result['mpi2']:>7.4f}s  {t_s/result['mpi2']:>7.2f}x")
    print(f"  {'MPI (4P)':<15} {result['mpi4']:>7.4f}s  {t_s/result['mpi4']:>7.2f}x")
    print(f"  {'MPI (8P)':<15} {result['mpi8']:>7.4f}s  {t_s/result['mpi8']:>7.2f}x")

