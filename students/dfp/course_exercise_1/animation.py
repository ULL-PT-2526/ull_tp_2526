import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.animation import FuncAnimation

plt.rcParams.update({
    'axes.labelcolor': 'white',
    'xtick.color': 'white',
    'ytick.color': 'white',
    'axes.titlecolor': 'white',
    'legend.facecolor': 'black',
    'legend.edgecolor': 'white',
    'legend.fontsize': 'medium',
    'text.color': 'white',
    'figure.facecolor': 'black',
    'figure.edgecolor': 'white',
    'axes.facecolor': 'black',
    'axes.edgecolor': 'white'
})

suffix = '_disk'
# === Configuration ===
filename = f"output{suffix}.dat"
interval_ms = 50   # delay between frames
point_size = 10    # particle marker size


# === Load the data ===
data = np.loadtxt(filename)

# Determine number of particles
ncols = data.shape[1]
n_particles = (ncols - 1) // 3
times = data[:, 0]

# Extract positions: shape (frames, n_particles, 3)
positions = data[:, 1:].reshape(len(times), n_particles, 3)

# === Setup 3D plot ===
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection='3d')
ax.set_title("N-body Simulation")

scat = ax.scatter([], [], [], s=point_size, color='w', alpha=0.75)

pos_min = positions.min()
pos_max = positions.max()
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_zlim(-3, 3)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.xaxis.pane.set_facecolor((0.1, 0.1, 0.1, 0.8))
ax.yaxis.pane.set_facecolor((0.1, 0.1, 0.1, 0.8))
ax.zaxis.pane.set_facecolor((0.1, 0.1, 0.1, 0.8))

# === Animation update function ===
def update(frame):
    scat._offsets3d = (
        positions[frame, :, 0],
        positions[frame, :, 1],
        positions[frame, :, 2],
    )
    ax.set_title(f"t = {times[frame]:.2f}")
    return scat,


# === Run animation ===
ani = FuncAnimation(fig, update, frames=len(times), interval=interval_ms, blit=False)

# Save animation as GIF
gif_filename = f"animation{suffix}.gif"
ani.save(gif_filename, writer='pillow', fps=1000/interval_ms)
print(f"Animation saved as {gif_filename}")

plt.show()
