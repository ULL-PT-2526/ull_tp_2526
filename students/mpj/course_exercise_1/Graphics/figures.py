from matplotlib.animation import FuncAnimation
import matplotlib.pylab as plt
import pandas as pd

# Reading the output file
path = "C:/Users/uSer/Documents/Máster Astrofísica/Segundo curso/Programación/fortran_course/students/mpj/course_exercise_1/"
df = pd.read_csv(path + 'output.dat', sep=r"\s+", header=None)  # separador por espacios
df.columns = ["x", "y", "z"]

n_part = 3
steps = int(len(df) / 3)


# =========================================================
# PARTICLES MOVEMENT ANIMATION
# =========================================================

# Colors for each particle
colors = ['magenta', 'gold', 'silver']

# Initialising lists to store the positions of the particles
x_data = [[] for _ in range(n_part)]
y_data = [[] for _ in range(n_part)]

# Initialising the figure
fig, ax = plt.subplots(1, 1)
lines_trace = [ax.plot([], [], '-', color=colors[i])[0] for i in range(n_part)]
markers = [ax.plot([], [], 'o', color=colors[i])[0] for i in range(n_part)]

x_min, x_max = df["x"].min(), df["x"].max()
y_min, y_max = df["y"].min(), df["y"].max()

# 10% margin of the larger axis range to keep all particles visible
margin = 0.1 * max(x_max - x_min, y_max - y_min)

ax.set_xlim(x_min - margin, x_max + margin)
ax.set_ylim(y_min - margin, y_max + margin)
ax.set_aspect('equal', 'box')

ax.set_xlabel("x")
ax.set_ylabel("y")

title = ax.set_title('')


# Initialisation function
def init():
    for line, marker in zip(lines_trace, markers):
        line.set_data([], [])
        marker.set_data([], [])
    return lines_trace + markers


# function that creates each frame of the animation
def update(frame):
    for i in range(n_part):
        idx = frame*n_part + i
        x, y = df.loc[idx, ["x", "y"]]
        x_data[i].append(x)
        y_data[i].append(y)
        markers[i].set_data([x], [y])
        lines_trace[i].set_data(x_data[i], y_data[i])

    title.set_text(f"t = {frame}")

    return lines_trace, markers, title


# Creating the animation
ani = FuncAnimation(fig, update, init_func=init, frames=steps, interval=60, blit=False)
#ani.save(path+"Graphics/animation.gif", writer="pillow", fps=10)

plt.show()
