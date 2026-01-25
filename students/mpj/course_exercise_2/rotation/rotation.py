import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ==========================
# Parámetros
# ==========================
filename = "/home/josemiguelmp/fortran_course/students/mpj/course_exercise_2/output_mpi.dat"
gif_name = "/home/josemiguelmp/fortran_course/students/mpj/course_exercise_2/animations/rotation_mpi.gif"
limit = 8       # Zoom de la cámara
skip = 2        # Dibuja uno de cada dos frames

# ==========================
# Leer datos
# ==========================
data = np.loadtxt(filename)

data = data[::skip, :]   # Aplicamos el salto de frames

time = data[:, 0]
coords = data[:, 1:]

# Calculamos N basándonos en que cada partícula tiene (x, y, z)
n_particles = (data.shape[1] - 1) // 3
positions = data[:, 1:].reshape(len(time), n_particles, 3)

# ==========================
# Figura 3D
# ==========================
fig = plt.figure(figsize=(10, 10), facecolor='black')
ax = fig.add_subplot(projection='3d', facecolor='black')

# Central body
scat_bh = ax.scatter([], [], [], s=80, c='white', edgecolors='red', lw=1, zorder=10)

scat_stars = ax.scatter([], [], [], s=2, c='cyan', alpha=0.4, edgecolors='none')

# Estética del espacio
ax.grid(False)
ax.set_axis_off() 

def init():
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    return scat_stars, scat_bh

def update(frame):
    # Posiciones del núcleo (BH)
    bx, by, bz = positions[frame, 0, :]
    scat_bh._offsets3d = ([bx], [by], [bz])
    
    # Posiciones de las estrellas
    sx = positions[frame, 1:, 0]
    sy = positions[frame, 1:, 1]
    sz = positions[frame, 1:, 2]
    scat_stars._offsets3d = (sx, sy, sz)
    
    ax.set_title(f"Simulación Galáctica | t = {time[frame]:.2f}", 
                 color='white', fontsize=12, pad=-20)
    
    # Efecto de rotación de cámara lento para que mole más
    ax.view_init(elev=30, azim=frame * 0.5)
    
    return scat_stars, scat_bh

# ==========================
# Renderizado
# ==========================
ani = FuncAnimation(fig, update, frames=len(time), init_func=init, 
                    interval=40, blit=False)

print(f"Generando {gif_name} (esto puede tardar un poco)...")
ani.save(gif_name, writer="pillow", fps=25, dpi=80)
print(f"¡Hecho! Mira el archivo {gif_name}")

plt.show()