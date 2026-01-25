# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 11:41:38 2025

@author: angeles
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


data = np.loadtxt('output.dat')
t = data[:, 0]
n_particles = (data.shape[1] - 1) // 3  # We substract the time collumn and then we divide the data by 3 for the particles positions and detect
# different particles


# We separate the coordinates

positions = []
for i in range(n_particles):
    x = data[:, 1 + 3*i]
    y = data[:, 2 + 3*i]
    z = data[:, 3 + 3*i]
    positions.append((x, y, z))

fig, ax = plt.subplots()

# We create lines and points for each particle detected

lines = [ax.plot([], [], '-', color='grey', alpha=0.3)[0] for _ in range(n_particles)]
points = [ax.plot([], [], 'o')[0] for _ in range(n_particles)]


ax.set_xlabel("x")
ax.set_ylabel("y")

# axis limits
a_x = np.concatenate([p[0] for p in positions])
a_y = np.concatenate([p[1] for p in positions])
ax.set_xlim(a_x.min() - 1, a_x.max() + 1)
ax.set_ylim(a_y.min() - 1, a_y.max() + 1)

def animate(i):

    for k in range(n_particles):
        x, y, z = positions[k]
        lines[k].set_data(x[:i], y[:i])
        points[k].set_data(x[i], y[i])
    
    ax.set_title(r"$t = %.3f$" % t[i])
  
    
    return lines,points

    


ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=30)
    
ani.save('particles.gif',  writer='pillow')