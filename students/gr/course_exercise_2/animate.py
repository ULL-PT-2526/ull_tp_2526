import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

# --------- Paths ---------
n_particles = 100
file_characteristic = f"mpi_n{n_particles}"

script_dir = os.path.dirname(os.path.abspath(__file__))
anim_dir = os.path.join(script_dir, "animations")
os.makedirs(anim_dir, exist_ok=True)

filename = os.path.join(script_dir, "outputs", f"output_{file_characteristic}.dat")

output_file = os.path.join(anim_dir, f"animation_{file_characteristic}.gif")

# --------- Parameters ---------
interval = 150  # ms between frames
tail = 50  # trail length (0 = no trail)

# --------- Data loading ---------
data = np.loadtxt(filename)

time = data[:, 0]
pos = data[:, 1:]  # remove time column

# Number of particles
n_particles = pos.shape[1] // 3

# Reshape to (n_frames, n_particles, 3)
pos = pos.reshape(len(time), n_particles, 3)
n_frames = pos.shape[0]

# --------- Figure setup ---------
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# --------- Automatic limits from data ---------
xmin, ymin, zmin = pos.min(axis=(0, 1))
xmax, ymax, zmax = pos.max(axis=(0, 1))

"""pos0 = pos[0]

xmin, ymin, zmin = pos0.min(axis=0)
xmax, ymax, zmax = pos0.max(axis=0)"""

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)

# Equal aspect ratio
ax.set_box_aspect([1, 1, 1])

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# --------- Artists ---------
particle_color = "black"
particle_size = 2

points = [
    ax.plot(
        [],
        [],
        [],
        marker="o",
        linestyle="None",
        color=particle_color,
        markersize=particle_size,
    )[0]
    for _ in range(n_particles)
]

trails = [
    ax.plot([], [], [], linestyle="-", color=particle_color, alpha=0.3, linewidth=0.8)[
        0
    ]
    for _ in range(n_particles)
]


# --------- Animation functions ---------
def init():
    for p, t in zip(points, trails):
        p.set_data([], [])
        p.set_3d_properties([])
        t.set_data([], [])
        t.set_3d_properties([])
    return points + trails


def update(frame):
    for i in range(n_particles):
        x, y, z = pos[frame, i]
        points[i].set_data([x], [y])
        points[i].set_3d_properties([z])

        if tail > 0:
            start = max(0, frame - tail)
            xs = pos[start : frame + 1, i, 0]
            ys = pos[start : frame + 1, i, 1]
            zs = pos[start : frame + 1, i, 2]
            trails[i].set_data(xs, ys)
            trails[i].set_3d_properties(zs)

    ax.set_title(f"t = {time[frame]:.3f}")
    return points + trails


# --------- Create animation ---------
ani = FuncAnimation(
    fig, update, frames=n_frames, init_func=init, interval=interval, blit=False
)

# --------- Save animation ---------
ani.save(output_file, writer="pillow", dpi=150)
plt.close(fig)

print(f"Animation saved to: {output_file}")
