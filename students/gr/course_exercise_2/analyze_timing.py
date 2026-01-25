import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os


# Create folder to save the figures
script_dir = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(script_dir, "figures")
os.makedirs(FIG_DIR, exist_ok=True)


# ============================================================
# Read file
# ============================================================
def read_timing(filename="timing.dat"):

    df = pd.read_csv(
        filename, delim_whitespace=True, names=["mode", "N", "cores", "time"]
    )

    return df


# ============================================================
# Compute speedup and efficiency
# ============================================================
def compute_metrics(df):
    results = []

    for N in sorted(df["N"].unique()):
        dfN = df[df["N"] == N]

        # serial reference number
        serial = dfN[dfN["mode"] == "serial"]
        if len(serial) == 0:
            print(f"No hay referencia serial para N={N}")
            continue

        T_serial = serial["time"].mean()

        for _, row in dfN.iterrows():
            mode = row["mode"]
            P = row["cores"]
            T = row["time"]

            speedup = T_serial / T
            efficiency = speedup / P

            results.append(
                {
                    "N": N,
                    "mode": mode,
                    "cores": P,
                    "time": T,
                    "speedup": speedup,
                    "efficiency": efficiency,
                }
            )

    return pd.DataFrame(results)


# ============================================================
# Plots
# ============================================================
def plot_times(df):
    plt.figure(figsize=(6, 4))

    for mode in df["mode"].unique():
        d = df[df["mode"] == mode]
        plt.plot(d["N"], d["time"], "o-", label=mode)

    plt.xlabel("Number of particles N")
    plt.ylabel("Execution time (s)")
    plt.title("Execution time vs N")
    plt.xscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "time_vs_N.png"), dpi=150)
    plt.close()


def plot_speedup(df):
    plt.figure(figsize=(6, 4))

    for mode in ["openmp", "mpi"]:
        d = df[df["mode"] == mode]
        plt.plot(d["N"], d["speedup"], "o-", label=mode)

    plt.xlabel("Number of particles N")
    plt.ylabel("Speedup")
    plt.title("Parallel speedup")
    plt.xscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "speedup.png"), dpi=150)
    plt.close()


def plot_efficiency(df):
    plt.figure(figsize=(6, 4))

    for mode in ["openmp", "mpi"]:
        d = df[df["mode"] == mode]
        plt.plot(d["N"], d["efficiency"], "o-", label=mode)

    plt.xlabel("Number of particles N")
    plt.ylabel("Efficiency")
    plt.title("Parallel efficiency")
    plt.ylim(0, 1.05)
    plt.xscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "efficiency.png"), dpi=150)
    plt.close()


# ============================================================
# PRINT SUMMARY
# ============================================================
def print_summary(df):
    print("\n==== SUMMARY ====\n")
    print(df.sort_values(["N", "mode", "cores"]))

    print("\nMean for mode:")
    print(df.groupby("mode")[["time", "speedup", "efficiency"]].mean())


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    filename = os.path.join(script_dir, "timing.dat")
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    df = read_timing(filename)

    metrics = compute_metrics(df)

    print_summary(metrics)

    plot_times(metrics)
    plot_speedup(metrics)
    plot_efficiency(metrics)

    # metrics.to_csv(os.path.join(FIG_DIR, "timing_analysis.csv"), index=False)

    print("\nGenerated files are in the 'figures' folder:")
    print(" - time_vs_N.png")
    print(" - speedup.png")
    print(" - efficiency.png")
