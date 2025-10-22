import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  

data = np.loadtxt("output.dat")

# Number of particles
n_particles = (data.shape[1] - 1) // 3 

# Extract times and positions
t = data[:, 0]
positions = data[:, 1:].reshape(len(t), n_particles, 3)

# Set up the 3D figure
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

# Define colors for particles and their trails
# It gives each particle a color from a list, depending on the total number of them
colors = plt.cm.tab10(np.linspace(0, 1, n_particles))

scat = ax.scatter(
    positions[0, :, 0],
    positions[0, :, 1],
    positions[0, :, 2],
    s=50,
    depthshade=False,
    c=colors,
)


x_min, x_max = positions[:, :, 0].min(), positions[:, :, 0].max()
y_min, y_max = positions[:, :, 1].min(), positions[:, :, 1].max()
z_min, z_max = positions[:, :, 2].min(), positions[:, :, 2].max()
ax.set_xlim(x_min - 0.5, x_max + 0.5)
ax.set_ylim(y_min - 0.5, y_max + 0.5)
ax.set_zlim(z_min - 0.5, z_max + 0.5)

# Trail length
trail_length = 25

# We create empty line objects for trails 
trails = [ax.plot([], [], [], lw=1, color=colors[i])[0] for i in range(n_particles)]


def update(frame):
    # Update 3D scatter positions
    scat._offsets3d = (
        positions[frame, :, 0],
        positions[frame, :, 1],
        positions[frame, :, 2],
    )
    ax.set_title(f"t = {t[frame]:.2f}")

    # Update trails for each particle
    for i in range(n_particles):
        start_frame = max(0, frame - trail_length)
        trail_data = positions[start_frame : frame + 1, i]
        trails[i].set_data(trail_data[:, 0], trail_data[:, 1])
        trails[i].set_3d_properties(trail_data[:, 2])

    return ()


ani = FuncAnimation(fig, update, frames=len(t), interval=50, blit=False)

# Sometimes we just want to see results on screen, not saving it...
save_choice = input("Do you want to save the animation as a video? [y/n]: ").strip().lower()

if save_choice == "y":
    print("Saving animation to animation.mp4... (this may take a while)")
    ani.save("animation.mp4", writer="ffmpeg", dpi=200)
    print("Video successfully saved in the current folder.")

plt.tight_layout()
plt.show()









