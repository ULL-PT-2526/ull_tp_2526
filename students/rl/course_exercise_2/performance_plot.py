import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv(
    "timings.txt",
    sep=r"\s+",
    comment="#",
    names=["N", "version", "workers", "time"]
)

# Select configurations
serial = df[(df["version"] == "serial") & (df["workers"] == 1)]

openmp_2 = df[(df["version"] == "openmp") & (df["workers"] == 2)]
openmp_4 = df[(df["version"] == "openmp") & (df["workers"] == 4)]

mpi_2 = df[(df["version"] == "mpi") & (df["workers"] == 2)]
mpi_4 = df[(df["version"] == "mpi") & (df["workers"] == 4)]

# Sort by N 
serial = serial.sort_values("N")
openmp_2 = openmp_2.sort_values("N")
openmp_4 = openmp_4.sort_values("N")
mpi_2 = mpi_2.sort_values("N")
mpi_4 = mpi_4.sort_values("N")


# Plot
plt.figure(figsize=(8, 6))

plt.plot(serial["N"], serial["time"], marker="o", label="Serial")

plt.plot(openmp_2["N"], openmp_2["time"], marker="s", linestyle="--",
         label="OpenMP (2 threads)")
plt.plot(openmp_4["N"], openmp_4["time"], marker="s",
         label="OpenMP (4 threads)")

plt.plot(mpi_2["N"], mpi_2["time"], marker="^", linestyle="--",
         label="MPI (2 processes)")
plt.plot(mpi_4["N"], mpi_4["time"], marker="^",
         label="MPI (4 processes)")

plt.xlabel("Number of particles (N)")
plt.ylabel("Execution time [s]")
plt.legend()

plt.xscale("log")
plt.yscale("log")

plt.tight_layout()
plt.savefig("time_plot.png", dpi=300)
plt.show()