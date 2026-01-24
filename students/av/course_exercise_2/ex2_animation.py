import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import argparse

import warnings
warnings.filterwarnings("ignore")



#--- Argument configuration ---#
parser = argparse.ArgumentParser(description="N-body simulation animation from an input file")

# Defining flag -i for an input file
parser.add_argument('-i', '--input', dest='input_file', type=str, default='output.dat',
                    help="Input data file (default: output.dat)")
parser.add_argument('-o', '--output', dest='output_file', type=str, default='ex1_animation',
                    help="Output animation name (Defatult: ex1_animation)")

args = parser.parse_args()
input_filename = args.input_file
output_filename = args.output_file



#--- Loading data ---#
print(f"Loading data from {input_filename}...")

# Reading the file
try:
    data = np.loadtxt(input_filename)
except FileNotFoundError:
    print("Simulation not created due to an error:")
    print(f"{input_filename} file not found.")
    exit()

# Obtaining time form first column
time = data[:, 0]

# Following columns are particle positions
positions_flat = data[:, 1:]

# Determining time steps
num_timesteps = data.shape[0]

# Managing variable particle number
total_columns = data.shape[1]

if (total_columns - 1) % 3 != 0:
    print(f"Error: '{input_filename}' file has unexpected format.")
    print(f"Waiting for 1 + (N*3) columnn but given just {total_columns}.")
    exit()

# Three columns (x,y,z) for each particle
# This give us the particle number
num_particles = (total_columns - 1) // 3

# Reshaping the array
positions = positions_flat.reshape(num_timesteps, num_particles, 3)

print(f"Data loaded: {num_timesteps} timesteps for {num_particles} particles detected.")



#--- 3D animation visualization configuration ---#

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Fixed axis limits from data
x_min, x_max = np.min(positions[:, :, 0]), np.max(positions[:, :, 0])
y_min, y_max = np.min(positions[:, :, 1]), np.max(positions[:, :, 1])
z_min, z_max = np.min(positions[:, :, 2]), np.max(positions[:, :, 2])

ax.set_xlim(x_min * 1.2, x_max * 1.2)
ax.set_ylim(y_min * 1.2, y_max * 1.2)
ax.set_zlim(z_min * 1.2, z_max * 1.2)

# Adding axis labels and title
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title(f"{num_particles} body simulation")

# Using a color map for particles (scalable)
colors = plt.cm.jet(np.linspace(0, 1, num_particles))

# For storing lines and dots
lines = []
points = []

# Loop for each particle (may be slow for lots of particles)
for i in range(num_particles):
    line_artist, = ax.plot([], [], [], '-', lw=1, color=colors[i], alpha=0.75)
    lines.append(line_artist)
    
    point_artist, = ax.plot([], [], [], 'o', markersize=6, color=colors[i], label=None)
    points.append(point_artist)

# Showing animation time
time_text = ax.text2D(0.05, 0.95, '', transform=ax.transAxes)

# Managing each particle path line length
TRAIL_LENGTH = 15 



#--- Needed functions to animate ---#

# Initializing the elements
def init():
    for line, point in zip(lines, points):
        line.set_data_3d([], [], [])
        point.set_data_3d([], [], [])
    time_text.set_text('')
    return lines + points + [time_text]

# Updating the frames for each particle
def animate(frame):
    
    start_index = max(0, frame - TRAIL_LENGTH + 1)

    for i in range(num_particles):
        
        x_trail = positions[start_index:frame+1, i, 0]
        y_trail = positions[start_index:frame+1, i, 1]
        z_trail = positions[start_index:frame+1, i, 2]
        
        lines[i].set_data_3d(x_trail, y_trail, z_trail)
        points[i].set_data_3d(x_trail[-1:], y_trail[-1:], z_trail[-1:])
    
    time_text.set_text(f'Time = {time[frame]:.2f} s')
    
    return lines + points + [time_text]



#--- Executing the animation ---#

print("Creating animation...")

ani = FuncAnimation(fig, animate, frames=num_timesteps,
                    init_func=init, blit=True, interval=30)

print("Saving animation...")
ani.save(f'{output_filename}.mp4', writer='ffmpeg', dpi=150)

#print("Showing animation...")
#plt.show()

print("Animation finished")