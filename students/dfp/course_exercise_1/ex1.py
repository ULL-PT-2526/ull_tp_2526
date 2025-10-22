import numpy as np
from time import time
from tqdm import tqdm

# ===============================================================
# 1. Read the setup file
# ===============================================================
def read_setup(filename):
    with open(filename, "r") as f:
        lines = [l for l in f if l.strip() and not l.strip().startswith("!")]
    dt = float(lines[0].split()[0])
    dt_out = float(lines[1].split()[0])
    t_end = float(lines[2].split()[0])
    n = int(lines[3].split()[0])
    data = np.loadtxt(lines[4:], dtype=float)
    m = data[:, 0]
    pos = data[:, 1:4]
    vel = data[:, 4:7]
    return dt, dt_out, t_end, n, m, pos, vel

# ===============================================================
# 2. Compute accelerations
# ===============================================================
def compute_accelerations(pos, m, eps=0.2):
    n = len(m)
    a = np.zeros_like(pos)
    for i in range(n):
        rji = pos - pos[i]
        r2 = np.sum(rji**2, axis=1) + eps**2
        r2[i] = np.inf
        r3 = r2 * np.sqrt(r2)
        a[i] = np.sum((m[:, None] * rji) / r3[:, None], axis=0)
    return a

# ===============================================================
# 3. Run the simulation
# ===============================================================
def nbody_simulation(setup_file="setup_disk.dat"):
    dt, dt_out, t_end, n, m, pos, vel = read_setup(setup_file)
    nstep = int(t_end / dt)

    acc = compute_accelerations(pos, m)
    vel += 0.5 * dt * acc

    start = time()

    for _ in tqdm(range(nstep), desc=f"Simulating ({n} particles)", unit="step"):
        pos += vel * dt
        acc = compute_accelerations(pos, m)
        vel += dt * acc

    elapsed = time() - start
    print(f"\nSimulation done in {elapsed:.3f} s "
          f"({n} particles, {nstep} steps).")

nbody_simulation("setup_disk.dat")
