import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D
import sys
import os

# --- visual configuration ---
subsample_step = 2    # skip frames to speed up playback
trail_length = 15     # trail length
fps = 30              # frames per second
size_central = 100    # size of the central black hole
size_particle = 10   # size of disk particles


# read command-line argument
if len(sys.argv) < 2:
    print("Usage: python3 plot_results.py <num_particles_input>")
    sys.exit(1)

n_input = int(sys.argv[1])

filename = f"output_disk_{n_input}.dat"

# read file header
with open(filename, 'r') as f:
    dt_val = float(f.readline())
    dt_out_val = float(f.readline())
    t_end_val = float(f.readline())
    disk_radius_val = float(f.readline())
    n_particles_total = int(f.readline())

print(f"--- Data: N={n_particles_total - 1}, Radius={disk_radius_val:.2f} ---")

# load numerical data
print("Loading numerical data...")
data = np.loadtxt(filename, skiprows=5)

# check dimensions
n_cols = data.shape[1]
n_calc = (n_cols - 1) // 3
if n_calc != n_particles_total:
    n_particles_total = n_calc

# extract and subsample
t_full = data[:, 0]
positions_full = data[:, 1:].reshape(len(t_full), n_particles_total, 3)

if subsample_step > 1:
    t = t_full[::subsample_step]
    positions = positions_full[::subsample_step, :, :]
else:
    t = t_full
    positions = positions_full

# figure setup
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title(f"Barnes-Hut (N={n_particles_total - 1})")

ax.grid(False)
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')

# colors for orbital particles (excluding index 0)
colors = plt.cm.viridis(np.linspace(0, 1, n_particles_total - 1))

# --- graphical objects ---
# central particle (index 0): black and large
scat_central = ax.scatter(
    positions[0, 0, 0],
    positions[0, 0, 1],
    positions[0, 0, 2],
    s=size_central,
    c='black',
    marker='o',
    depthshade=True,
    alpha=1.0,
    label='Central Mass'
)

# disk particles (indices 1 to end): colored and small
scat_orbit = ax.scatter(
    positions[0, 1:, 0],
    positions[0, 1:, 1],
    positions[0, 1:, 2],
    s=size_particle,
    c=colors,
    depthshade=False,
    alpha=0.7
)

# dynamic limits
limit_xy = disk_radius_val * 1.1
limit_z = disk_radius_val * 0.5
ax.set_xlim(-limit_xy, limit_xy)
ax.set_ylim(-limit_xy, limit_xy)
ax.set_zlim(-limit_z, limit_z)

# time label
time_text = ax.text2D(0.05, 0.95, '', transform=ax.transAxes, fontsize=12)

# trails 
show_trails = True
if n_particles_total > 800:
    show_trails = False

trails = []
if show_trails:
    trails = [
        ax.plot([], [], [], lw=0.8, color=colors[i], alpha=0.4)[0]
        for i in range(n_particles_total - 1)
    ]

def update(frame):
    # update central particle
    scat_central._offsets3d = (
        [positions[frame, 0, 0]],
        [positions[frame, 0, 1]],
        [positions[frame, 0, 2]]
    )

    # update disk particles (slice 1:)
    scat_orbit._offsets3d = (
        positions[frame, 1:, 0],
        positions[frame, 1:, 1],
        positions[frame, 1:, 2]
    )

    time_text.set_text(f"t = {t[frame]:.3f}")

    # update trails
    if show_trails:
        for i in range(n_particles_total - 1):
            idx_part = i + 1
            start = max(0, frame - trail_length)

            trail_data = positions[start:frame + 1, idx_part]
            trails[i].set_data(trail_data[:, 0], trail_data[:, 1])
            trails[i].set_3d_properties(trail_data[:, 2])

    if show_trails:
        return (scat_central, scat_orbit, time_text) + tuple(trails)
    else:
        return (scat_central, scat_orbit, time_text)

# save or display animation
print("Generating animation...")
ani = FuncAnimation(fig, update, frames=len(t), interval=1000 / fps, blit=False)

output_gif = f"animation_disk_{n_particles_total - 1}.gif"

save_choice = input(f"Save as '{output_gif}'? [y/n]: ").strip().lower()

if save_choice == "y":
    print("Saving GIF...")
    writer = PillowWriter(fps=fps)
    ani.save(output_gif, writer=writer)
    print(f"Done! Saved to: {output_gif}")
else:
    plt.tight_layout()
    plt.show()





