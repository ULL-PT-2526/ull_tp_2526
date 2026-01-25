import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ==========================
# Parameters
# ==========================
filename = "/home/josemiguelmp/fortran_course/students/mpj/course_exercise_2/output_mpi.dat"
gif_name = "/home/josemiguelmp/fortran_course/students/mpj/course_exercise_2/animations/spherical_collapse_mpi.gif"
limit = 1.2      # Initial radius of the sphere
skip = 1         # Change to 2 or more if the file is very heavy

# ==========================
# Reading data
# ==========================
print(f"Leyendo {filename}...")
data = np.loadtxt(filename)

if skip > 1:
    data = data[::skip, :]

time = data[:, 0]

# Calculating number of particles
# Each particle has coordinates (x, y, z)
n_particles = (data.shape[1] - 1) // 3
positions = data[:, 1:].reshape(len(time), n_particles, 3)

# ==========================
# 3D figure
# ==========================
fig = plt.figure(figsize=(8, 8), facecolor='black')
ax = fig.add_subplot(projection='3d', facecolor='black')
scat = ax.scatter([], [], [], s=3, c='yellow', alpha=0.6, edgecolors='none')

ax.grid(False)
ax.set_axis_off() 

def init():
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    return (scat,)

def update(frame):
    x = positions[frame, :, 0]
    y = positions[frame, :, 1]
    z = positions[frame, :, 2]
    
    scat._offsets3d = (x, y, z)
    ax.set_title(f"Colapso Esférico | t = {time[frame]:.3f}", 
                 color='white', fontsize=14, pad=-20)
    
    return (scat,)

# ==========================
# Render
# ==========================
ani = FuncAnimation(fig, update, frames=len(time), init_func=init, 
                    interval=30, blit=False)

print(f"Generando {gif_name}...")

ani.save(gif_name, writer="pillow", fps=30, dpi=100)
print(f"¡Hecho! Animación guardada en {gif_name}")

plt.show()