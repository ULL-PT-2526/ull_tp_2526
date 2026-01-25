import os
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

def generate_animation(filename, gif_filename, interval_ms=20, point_size=10, frame_step=1):
    try:
        data = np.loadtxt(filename)
    except Exception as e:
        print(f"Skipping {filename}: {e}")
        return False

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
    # Initialize trail lines for each particle
    trails = [ax.plot([], [], [], color='white', alpha=0.1, linewidth=1)[0] for _ in range(n_particles)]

    # Fixed bounds to keep view consistent across frames and runs
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_zlim(-3, 3)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.xaxis.pane.set_facecolor((0.1, 0.1, 0.1, 0.8))
    ax.yaxis.pane.set_facecolor((0.1, 0.1, 0.1, 0.8))
    ax.zaxis.pane.set_facecolor((0.1, 0.1, 0.1, 0.8))

    def update(frame):
        # Update particle positions
        scat._offsets3d = (
            positions[frame, :, 0],
            positions[frame, :, 1],
            positions[frame, :, 2],
        )

        # Update trails for each particle
        for i in range(n_particles):
            # Show trail from start up to current frame
            trail_x = positions[:frame+1, i, 0]
            trail_y = positions[:frame+1, i, 1]
            trail_z = positions[:frame+1, i, 2]
            trails[i].set_data_3d(trail_x, trail_y, trail_z)

        ax.set_title(f"t = {times[frame]:.2f}")
        return [scat] + trails

    frames = range(0, len(times), frame_step)
    ani = FuncAnimation(fig, update, frames=len(times), interval=interval_ms, blit=False)

    try:
        ani.save(gif_filename, writer='pillow', fps=1000/interval_ms)
        print(f"Animation saved as {gif_filename}")
        plt.close(fig)
        return True
    except Exception as e:
        print(f"Failed to save {gif_filename}: {e}")
        plt.close(fig)
        return False


if __name__ == "__main__":
    setups = ['hannu', 'disk', 'random']
    modes = ['mpi']

    for setup in setups: 
        for mode in modes:
            base = f"outputs/output_{setup}_{mode}.dat"
            filename = base

            if not os.path.exists(filename):
                print(f"Missing {base}, skipping.")
                continue

            gif_filename = f"animations/animation_{setup}_{mode}.gif"
            generate_animation(filename, gif_filename)

    # plt.show()
