import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os

# === FILE READING ===
filename = os.path.join(os.path.dirname(__file__), "output.dat")
data = np.loadtxt(filename)

time = data[:, 0]
nsteps, ncols = data.shape
nparticles = (ncols - 1) // 3  # calculate the number of particles in output.dat

# Reshape positions: shape (nsteps, nparticles, 3)
positions = data[:, 1:].reshape(nsteps, nparticles, 3)

# === COMPUTE GLOBAL LIMITS ===
x_min, y_min = positions[:, :, 0].min(), positions[:, :, 1].min()
x_max, y_max = positions[:, :, 0].max(), positions[:, :, 1].max()

# Add a 10% margin so particles stay within frame
margin_x = 0.1 * (x_max - x_min)
margin_y = 0.1 * (y_max - y_min)

# === FIGURE CONFIGURATION ===
fig, ax = plt.subplots(figsize=(7, 6))

colors = ["tab:red", "tab:blue", "tab:green", "tab:orange", "tab:purple", "tab:brown"]

# Initialize scatter points and trails
scatters = [
    ax.plot([], [], "o", color=colors[i % len(colors)], markersize=6)[0]
    for i in range(nparticles)
]
trails = [
    ax.plot([], [], "-", color=colors[i % len(colors)], linewidth=1, alpha=0.7)[0]
    for i in range(nparticles)
]

# === AXIS LIMITS (AUTO ZOOM) ===
ax.set_xlim(x_min - margin_x, x_max + margin_x)
ax.set_ylim(y_min - margin_y, y_max + margin_y)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Particle evolution over time")
ax.set_aspect("equal", adjustable="box")


# === ANIMATION FUNCTION ===
def update(frame):
    for i in range(nparticles):
        x, y = positions[: frame + 1, i, 0], positions[: frame + 1, i, 1]
        scatters[i].set_data(x[-1:], y[-1:])
        trails[i].set_data(x, y)
    ax.set_title(f"t = {time[frame]:.3f}")
    return scatters + trails


# === CREATE ANIMATION ===
ani = FuncAnimation(fig, update, frames=nsteps, interval=80, blit=False)

# === SAVE AS GIF ===
output_gif = os.path.join(os.path.dirname(__file__), "animation2D.gif")
writer = PillowWriter(fps=12)
ani.save(output_gif, writer=writer)
print(f"Animation saved.")

plt.show()
