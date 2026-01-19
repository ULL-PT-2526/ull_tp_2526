import numpy as np
import os

# ============================================================
# Global simulation parameters
# ============================================================
dt = 0.01
dt_out = 0.1
t_end = 2.0

# Initial sphere radius
R = 1.0

# Gravitational constant (set to 1 for convenience)
G = 1.0

# Particle mass
mass_per_particle = 1.0

# Cases to generate
Ns = [10, 100, 1000, 5000, 10000]

output_dir = "./students/gr/course_exercise_2/setups"
os.makedirs(output_dir, exist_ok=True)


# ============================================================
# Functions
# ============================================================
def random_in_sphere(N, R):

    u = np.random.rand(N)
    costheta = 2 * np.random.rand(N) - 1
    phi = 2 * np.pi * np.random.rand(N)

    r = R * u ** (1 / 3)
    theta = np.arccos(costheta)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z


def random_velocities(N):

    vx = np.random.randn(N)
    vy = np.random.randn(N)
    vz = np.random.randn(N)

    # Remove any net bulk motion
    vx -= vx.mean()
    vy -= vy.mean()
    vz -= vz.mean()

    return vx, vy, vz


def virialize_velocities(vx, vy, vz, N, m, R):

    M = N * m

    # Potential energy of a uniform sphere
    U = -3.0 / 5.0 * G * M * M / R

    # Current kinetic energy
    K_current = 0.5 * m * np.sum(vx**2 + vy**2 + vz**2)

    # Target kinetic energy from virial theorem
    K_target = abs(U) / 2.0

    # Rescaling factor
    factor = np.sqrt(K_target / K_current)

    vx *= factor
    vy *= factor
    vz *= factor

    return vx, vy, vz


# ============================================================
# File generation
# ============================================================
for N in Ns:

    filename = os.path.join(output_dir, f"setup_{N}.dat")

    # Positions
    x, y, z = random_in_sphere(N, R)

    # Initial random velocities
    vx, vy, vz = random_velocities(N)

    # Virialize the system
    vx, vy, vz = virialize_velocities(vx, vy, vz, N, mass_per_particle, R)

    with open(filename, "w") as f:

        # Header
        f.write(f"{dt}\n")
        f.write(f"{dt_out}\n")
        f.write(f"{t_end}\n")
        f.write(f"{N}\n")

        # Particles
        for i in range(N):
            f.write(
                f"{mass_per_particle} "
                f"{x[i]} {y[i]} {z[i]} "
                f"{vx[i]} {vy[i]} {vz[i]}\n"
            )

    print(f"Generated {filename}")

print("\nAll setups created in folder:", output_dir)
