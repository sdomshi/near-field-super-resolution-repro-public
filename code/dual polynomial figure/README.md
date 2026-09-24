# Lifted dual-polynomial localization experiment

This directory reproduces the finite Jacobi--Anger lifted SDP experiment in
the paper. The code generates data with the quadratic Fresnel phase model,
solves the robust dual SDP, extracts the active range--angle support, refits
the complex amplitudes, and exports the paper-facing figure and diagnostics.

## Requirements

- MATLAB R2019b or later
- [CVX 2.2](http://cvxr.com/cvx/)
- SDPT3 4.0 configured through CVX

The experiment does not require external data files.

## Run

Start MATLAB, initialize CVX if necessary, change to this directory, and run:

```matlab
result = run_lifted_dual_figure1_repro();
```

The script is a function and can be called from any working directory. Output
paths are resolved relative to the function file, not the caller's directory.

Typical runtime with SDPT3 is approximately one minute, but it depends on the
MATLAB, CVX, BLAS, and solver versions. Compare numerical values using the
tolerances below rather than bitwise equality.
