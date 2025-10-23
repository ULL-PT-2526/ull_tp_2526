import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


# Names of the files
input = "output.dat" # ex1 output
output = "simulation.gif" # name of the animation output from this script

# Each line of the ex1 output has been formatted to have the form:
# time p1x p1y p1z p2x p2y p2z ... pnx pny pnz
data = np.loadtxt(input)
time = data[:, 0] # shape: [frames]
pos = data[:, 1:] # shape: [frames, 1+3*n], where n is the num of particles

n = pos.shape[1]//3 # number of particles

# Structure the positions a little better/clearer:
pos = pos.reshape(len(time), n, 3) # now pos[frame, particle, coordinate]
# example: pos[3, 0, 1] will be the x coordinate of particle 1 at frame 4


# With this, we can set up the 2D plot
fig, ax = plt.subplots(figsize=(6, 6))

ax.set_aspect('equal', 'box')
ax.set_xlabel("x")
ax.set_ylabel("y") 
ax.minorticks_on()

# determine the limits by flattening and then finding min. and max. values
x_flat = pos[:, :, 0].flatten()
y_flat = pos[:, :, 1].flatten()
margin = 0.1*max(x_flat.max() - x_flat.min(), y_flat.max() - y_flat.min())
ax.set_xlim(x_flat.min() - margin, x_flat + margin)
ax.set_ylim(y_flat.min() - margin, y_flat.max() + margin)

# set up elements of the plot
points, = ax.plot([], [], 'o', color='k', markersize=6)
paths = [ax.plot([], [], '-', lw=0.8, color='gray', alpha=0.4)[0] for _ in range(n)]
title = ax.text(0.5, 1.02, '', transform=ax.transAxes, fontsize=12, ha='center') # adjusted position and added horixontal alignment

def update(frame):
    # print how far along it is
    if frame % 100 == 0: # every 100 frames
        print("Frame {}/{}".format(frame, len(time)))
    
    # update particle positions
    x = pos[frame, :, 0]
    y = pos[frame, :, 1]
    points.set_data(x, y) # fill in what we set u earlier

    # update trajectories
    for i in range(n):
        paths[i].set_data(pos[:frame+1, i, 0], pos[:frame+1, i, 1]) # all x and y coords from the start until this frame

    # update time (in the title)
    title.set_text(f"t = {time[frame]:.2f}")
    
    return [points, *paths, title]

# Finally, create the animation
anim = FuncAnimation(fig, update, frames=len(time), interval=50, blit=True)

# and save it as a GIF
anim.save(output, writer=PillowWriter(fps=20))

plt.close()
print("Animation saved as {}".format(output))