import random
import sys

def create_nbody_input(n_particles, filename):
    # Simulation parameters
    dt = 0.01      # Time step 
    dt_out = 1.0   # Output frequency 
    t_end = 100.0    # Total time
    
    with open(filename, 'w') as f:
        # Write header
        f.write(f"{dt}\n")
        f.write(f"{dt_out}\n")
        f.write(f"{t_end}\n")
        f.write(f"{n_particles}\n")
        
        # Generate random particles
        for i in range(n_particles):
            mass = random.random()
            if i == 0: mass = 1.0
            # Random position in a 3D cube [-1, 1]
            px, py, pz = [random.uniform(-1, 1) for _ in range(3)]
            # Random small initial velocities
            vx, vy, vz = [random.uniform(-0.05, 0.05) for _ in range(3)]
            
            f.write(f"{mass:.6f} {px:.6f} {py:.6f} {pz:.6f} {vx:.6f} {vy:.6f} {vz:.6f}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_input.py <number_of_particles>")
    else:
        n = int(sys.argv[1])
        outfile = f"input_{n}.dat"
        create_nbody_input(n, outfile)
        print(f"Successfully created {outfile} with {n} particles.")