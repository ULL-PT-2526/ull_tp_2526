import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os

# === FILE READING ===
filename = os.path.join(os.path.dirname(__file__), "output.dat")
data = np.loadtxt(filename)

time = data[:, 0]
nsteps, ncols = data.shape
nparticles = (ncols - 1) // 3

# Reshape positions: shape (nsteps, nparticles, 3)
positions = data[:, 1:].reshape(nsteps, nparticles, 3)

# === COMPUTE GLOBAL LIMITS ===
xyz_min = positions.min(axis=(0, 1))
xyz_max = positions.max(axis=(0, 1))

x_min, y_min, z_min = xyz_min
x_max, y_max, z_max = xyz_max

# Add a 10% margin so particles stay within frame
margin_x = 0.1 * (x_max - x_min)
margin_y = 0.1 * (y_max - y_min)

# === FIGURE CONFIGURATION ===
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")

colors = ["tab:red", "tab:blue", "tab:green", "tab:orange", "tab:purple", "tab:brown"]

# Initialize scatter points and trails
scatters = [
    ax.plot([], [], [], "o", color=colors[i % len(colors)], markersize=6)[0]
    for i in range(nparticles)
]
trails = [
    ax.plot([], [], [], "-", color=colors[i % len(colors)], linewidth=1, alpha=0.7)[0]
    for i in range(nparticles)
]

# === AXIS LIMITS (AUTO ZOOM) ===
ax.set_xlim(x_min - margin_x, x_max + margin_x)
ax.set_ylim(y_min - margin_y, y_max + margin_y)
ax.set_zlim(-1, 1)  # practically flat z-axis

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Particle evolution over time")


# === ANIMATION FUNCTION ===
def update(frame):
    for i in range(nparticles):
        x, y, z = (
            positions[: frame + 1, i, 0],
            positions[: frame + 1, i, 1],
            positions[: frame + 1, i, 2],
        )
        scatters[i].set_data(x[-1:], y[-1:])
        scatters[i].set_3d_properties(z[-1:])
        trails[i].set_data(x, y)
        trails[i].set_3d_properties(z)
    ax.set_title(f"t = {time[frame]:.3f}")
    return scatters + trails


# === CREATE ANIMATION ===
ani = FuncAnimation(fig, update, frames=nsteps, interval=80, blit=False)

# === SAVE AS GIF ===
output_gif = os.path.join(os.path.dirname(__file__), "animation3D.gif")
writer = PillowWriter(fps=12)
ani.save(output_gif, writer=writer)
print(f"Animation saved.")

plt.show()
