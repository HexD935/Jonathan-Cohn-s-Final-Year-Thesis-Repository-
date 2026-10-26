# Wavelet-Coefficient Impedance Estimation (EEE4022S, 2026)

MATLAB/Simulink code for the thesis "Direct Wavelet-Coefficient Impedance
Estimation in Unbalanced Three-Phase Four-Wire Networks", J. P. Cohn,
University of Cape Town.

| Folder | Thesis section |
|---|---|
| `step05_dwt/` | Sym4 decomposition and reconstruction (Sec. 3.5) |
| `step06-08_load_impedance/` | Z_Load from wavelet coefficients (Secs. 3.6–3.8) |
| `step09_thevenin/` | Z_th from detected load switches (Sec. 3.9) |
| `step10_three_phase/` | Three-phase model and neutral impedance (Sec. 3.10) |

## How to run
1. Run the Simulink model for the step (10 s, fs = 8192 Hz). It writes `out.simout`.
2. Run the matching MATLAB script. It prints the error summary and plots the figures.
