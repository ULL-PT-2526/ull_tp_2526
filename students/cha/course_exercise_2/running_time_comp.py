import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Number of particles used in the simulations
N_values = np.array([50, 500, 2000, 10000])

# Execution times (in seconds) for Barnes–Hut implementations
times_serial   = np.array([0.63, 11.91, 57.47, 380.82]) 
times_omp_4    = np.array([0.75, 8.71, 33.69, 198.56]) 
times_omp_8    = np.array([4.68, 14.95, 41.22, 184.31]) 
times_mpi_1    = np.array([1.08, 12.41, 57.73, 350.42]) 
times_hybrid_1 = np.array([0.92, 11.03, 39.21, 245.23]) 
times_mpi_2    = np.array([0.98, 13.20, 54.40, 341.87]) 
times_hybrid_2 = np.array([0.85, 10.02, 38.49, 242.02]) 


# Execution times for the direct N-body method (serial O(N^2))
N_nbody = np.array([50, 500, 2000])
times_nbody = np.array([0.49, 40.78, 647.83])

# Barnes–Hut serial times corresponding to the same N values
times_bh_serial = times_serial[:3]


# Theoretical scaling models used for curve fitting
def model_n2(N, a):
    return a * N**2

def model_nlogn(N, b):
    return b * N * np.log2(N)


# Fit the theoretical models to the measured data using SciPy
(a_fit,), _ = curve_fit(model_n2, N_nbody, times_nbody)
(b_fit,), _ = curve_fit(model_nlogn, N_nbody, times_bh_serial)

# Smooth range of N values for plotting the fitted curves
N_fit = np.logspace(np.log10(50), np.log10(10000), 300)


# ============================================================
# Figure 1: Direct N-body vs Barnes–Hut (Serial) with fits
# ============================================================

plt.figure(figsize=(8, 6))

# Measured execution times
plt.plot(N_nbody, times_nbody, 'o', color='magenta',
         markersize=7, label='Direct N-body')

plt.plot(N_nbody, times_bh_serial, 's', color='black',
         markersize=7, label='Barnes-Hut (serial)')

# Fitted theoretical scaling laws
plt.plot(N_fit, model_n2(N_fit, a_fit), '--', color='magenta',
         linewidth=2.5, alpha=0.7, label=r'Fit $\sim N^2$')

plt.plot(N_fit, model_nlogn(N_fit, b_fit), '--', color='black',
         linewidth=2.5, alpha=0.7, label=r'Fit $\sim N\log N$')

# Log–log scaling highlights the asymptotic behavior
plt.xscale('log')
plt.yscale('log')

plt.xlabel('Number of particles (N)', fontsize=12)
plt.ylabel('Execution time (seconds)', fontsize=12)
plt.title('Direct N-body vs Barnes-Hut (Serial)', fontsize=14, fontweight='bold')

plt.legend()

plt.tight_layout()
plt.savefig('nbody_vs_bh_serial.png', dpi=300)
plt.show()


# ============================================================
# Figure 2: Barnes–Hut performance analysis (time & speedup)
# ============================================================

# Compute speedups relative to the serial Barnes–Hut version
speedup_omp_4    = times_serial / times_omp_4
speedup_omp_8    = times_serial / times_omp_8
speedup_mpi_1    = times_serial / times_mpi_1
speedup_hybrid_1 = times_serial / times_hybrid_1
speedup_mpi_2    = times_serial / times_mpi_2    
speedup_hybrid_2 = times_serial / times_hybrid_2 

# Create a figure with two subplots: execution time and speedup
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), sharex=True)

fig.suptitle('Barnes–Hut Performance Analysis',
             fontsize=18, fontweight='bold', y=0.98)

# Execution time comparison
ax1.plot(N_values, times_serial, 'o--', color='black',
         label='Serial', linewidth=2)
ax1.plot(N_values, times_omp_4,  's-',  color='red',
         label='OpenMP (4 threads)', linewidth=2)
ax1.plot(N_values, times_omp_8,  '^:',  color='orange',
         label='OpenMP (8 threads)', linewidth=1.5)

ax1.plot(N_values, times_mpi_1, 'D-', color='blue',
         label='MPI v1 (4 processes)', linewidth=1.5, alpha=0.6)
ax1.plot(N_values, times_mpi_2, 'D-', color='darkblue',
         label='MPI v2 (4 processes, optimized)', linewidth=2.5)

ax1.plot(N_values, times_hybrid_1, 'x-.', color='limegreen',
         label='Hybrid v1 (2p x 2t)', linewidth=1.5, alpha=0.6)
ax1.plot(N_values, times_hybrid_2, 'x-', color='darkgreen',
         label='Hybrid v2 (2p x 2t, optimized)', linewidth=2.5)

ax1.set_xscale('log')
ax1.set_xlabel('Number of particles (N)', fontsize=12)
ax1.set_ylabel('Execution time (seconds)', fontsize=12)
ax1.legend(fontsize=9, loc='upper left')

# Speedup comparison
ax2.plot(N_values, speedup_omp_4, 's-', color='red',
         label='Speedup OpenMP (4 threads)', linewidth=2)
ax2.plot(N_values, speedup_omp_8, '^:', color='orange',
         label='Speedup OpenMP (8 threads)', linewidth=1.5)

ax2.plot(N_values, speedup_mpi_1, 'D-', color='blue',
         label='Speedup MPI v1', linewidth=1.5, alpha=0.6)
ax2.plot(N_values, speedup_mpi_2, 'D-', color='darkblue',
         label='Speedup MPI v2', linewidth=2.5)

ax2.plot(N_values, speedup_hybrid_1, 'x-.', color='limegreen',
         label='Speedup Hybrid v1', linewidth=1.5, alpha=0.6)
ax2.plot(N_values, speedup_hybrid_2, 'x-', color='darkgreen',
         label='Speedup Hybrid v2', linewidth=2.5)

# Reference line: speedup = 1 (serial baseline)
ax2.axhline(y=1, color='black', linestyle='--', alpha=0.8)

ax2.set_xlabel('Number of particles (N)', fontsize=12)
ax2.set_ylabel('Speedup factor', fontsize=12)
ax2.legend(fontsize=9, loc='upper left')

# Highlight the maximum speedup achieved at the largest N
all_speedups = [
    speedup_omp_4[-1], speedup_omp_8[-1],
    speedup_mpi_1[-1], speedup_mpi_2[-1],
    speedup_hybrid_1[-1], speedup_hybrid_2[-1]
]
best_speedup = max(all_speedups)

ax2.annotate(f"Max speedup: {best_speedup:.2f}×",
             (N_values[-1], best_speedup),
             xytext=(-60, -20), textcoords='offset points',
             arrowprops=dict(arrowstyle="->"),
             fontweight='bold', color='orange')

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('performance_comparison.png', dpi=300)
plt.show()