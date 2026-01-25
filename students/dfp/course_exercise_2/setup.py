import math
import random
from typing import List, Tuple

# ===============================================================
# Writer for the setup file
# ===============================================================
def write_setup(filename: str, dt: float, dt_out: float, t_end: float,
                particles: List[Tuple[float, float, float, float, float, float, float]]):
    with open(filename, 'w') as f:
        f.write(f"{dt}         ! dt (time step)\n")
        f.write(f"{dt_out}      ! dt_out (output interval)\n")
        f.write(f"{t_end}       ! t_end (total simulation time)\n")
        f.write(f"{len(particles)}          ! n (number of particles)\n")
        for p in particles:
            m, x, y, z, vx, vy, vz = p
            f.write(f"{m:.16e} {x:.16e} {y:.16e} {z:.16e} {vx:.16e} {vy:.16e} {vz:.16e}\n")


# ===============================================================
# 1. Disk Galaxy
# ===============================================================
def disk_galaxy(n: int, radius: float = 2.0, thickness: float = 0.1, mass: float = 1.0,
                seed: int = None, eps: float = 0.3) -> List[Tuple[float, float, float, float, float, float, float]]:
    rnd = random.Random(seed)
    particles = []

    # Enclosed mass for exponential disk
    def enclosed_mass(R: float) -> float:
        return n * mass * (1 - (1 + R/radius) * math.exp(-R/radius))

    for _ in range(n):
        R = -radius * math.log(1 - rnd.random())
        theta = rnd.random() * 2 * math.pi
        x = R * math.cos(theta)
        y = R * math.sin(theta)
        z = rnd.gauss(0, thickness)

        # Circular velocity using enclosed mass
        M_enc = enclosed_mass(R)
        v_c = math.sqrt(M_enc / max(R + eps, 1e-6))

        # Small velocity noise (0.5% of circular speed)
        vt = v_c * (1 + rnd.gauss(0, 0.005))
        vx = -vt * math.sin(theta)
        vy = vt * math.cos(theta)
        vz = rnd.gauss(0, 0.005 * vt)

        particles.append((mass, x, y, z, vx, vy, vz))

    return particles


# ===============================================================
# 2. Random uniform box
# ===============================================================
def random_uniform(n: int, scale: float = 2.0, mass: float = 1.0, seed: int = None):
    rnd = random.Random(seed)
    particles = []
    for _ in range(n):
        x = rnd.uniform(-scale, scale)
        y = rnd.uniform(-scale, scale)
        z = rnd.uniform(-scale, scale)
        vx = rnd.uniform(-0.1, 0.1)
        vy = rnd.uniform(-0.1, 0.1)
        vz = rnd.uniform(-0.1, 0.1)
        particles.append((mass, x, y, z, vx, vy, vz))
    return particles


# ===============================================================
# Generator
# ===============================================================
def generate_particles(ptype: str, n: int, scale: float, mass: float, seed: int = None):
    if ptype == 'disk':
        return disk_galaxy(n, radius=scale, mass=mass, seed=seed)
    elif ptype == 'random':
        return random_uniform(n, scale=scale, mass=mass, seed=seed)


# ===============================================================
# Generate disk and random setups
# ===============================================================

n = 500
total_mass = 1.0
mass = total_mass / n
dt = 5e-5
dt_out = 0.05
t_end = 8.0
scale = 2.0
seed = 42

for t in ['disk', 'random']:
    particles = generate_particles(t, n, scale, mass, seed)
    write_setup(f"setup_{t}.dat", dt, dt_out, t_end, particles)
    print(f"Wrote {len(particles)} particles to setup_{t}.dat (type={t})")
