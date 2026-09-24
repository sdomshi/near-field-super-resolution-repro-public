# Reproducibility code for *A Mathematical Theory of Near-Field Super-Resolution*

**Author:** Sajad Daei  
**Contact:** [sajado@kth.se](mailto:sajado@kth.se), 
[sajaddaeiomshi@gmail.com](mailto:sajaddaeiomshi@gmail.com)

This repository contains the MATLAB code used for the numerical examples and
computational figures accompanying

> S. Daei, G. Fodor, and M. Skoglund,  
> *A Mathematical Theory of Near-Field Super-Resolution*.

The repository contains four primary, self-contained MATLAB functions. No
external input datasets are required.

## Repository contents

| Paper component | MATLAB entry point |
| --- | --- |
| Derivative-route QPAC support-class evaluation | [`derivative example/run_derivative_example_qpac.m`](derivative%20example/run_derivative_example_qpac.m) |
| Lag-correlation-support QPAC evaluation | [`lcs_example/run_lcs_example_qpac.m`](lcs_example/run_lcs_example_qpac.m) |
| Finite-harmonic lifted dual-polynomial experiment | [`dual_polynomial_figure/run_lifted_dual_figure1_repro.m`](dual_polynomial_figure/run_lifted_dual_figure1_repro.m) |
| Exact-kernel and QPAC-bound comparison plots | [`kernel bounds/run_bounds_true_kernel_repro.m`](kernel%20bounds/run_bounds_true_kernel_repro.m) |

Each experiment directory contains a separate README with its mathematical
scope, execution instructions, and generated outputs.

## Requirements

- MATLAB R2019b or later.
- The dual-polynomial experiment additionally requires:
  - [CVX 2.2](http://cvxr.com/cvx/);
  - SDPT3 4.0 configured as a CVX solver.
- The derivative-route, lag-correlation-support, and kernel-bound examples do
  not require CVX or external input files.

## Running the experiments

Start MATLAB in the repository root and record the root directory:

```matlab
repo_root = pwd;
```

### 1. Derivative-route QPAC example

```matlab
cd(fullfile(repo_root, 'derivative example'));
row = run_derivative_example_qpac();
```

The function evaluates the derivative-route QPAC quantities and writes the
paper-facing LaTeX tables and a MATLAB results file.

### 2. Lag-correlation-support QPAC example

```matlab
cd(fullfile(repo_root, 'lcs_example'));
result = run_lcs_example_qpac();
```

The function evaluates the lag-correlation-support construction and writes the
reported LaTeX table and MATLAB results file.

### 3. Finite-harmonic lifted dual-polynomial experiment

```matlab
cd(fullfile(repo_root, 'dual_polynomial_figure'));
result = run_lifted_dual_figure1_repro();
```

This function generates data from the quadratic Fresnel phase model, solves the
finite Jacobi--Anger lifted dual SDP, extracts the range--angle support, refits
the source amplitudes, and writes:

- `figures/dualpol1.pdf`, `dualpol1.png`, and `dualpol1.fig`;

### 4. Exact-kernel and QPAC-bound figures

```matlab
cd(fullfile(repo_root, 'kernel bounds'));
out = run_bounds_true_kernel_repro();
```

The function compares the exact channel magnitudes with the pointwise and
support-uniform bounds for the eight Hermite and derivative channels
`K`, `H`, `dK`, `dH`, `d2K`, `d3K`, `d2H`, and `d3H`.

## Interpretation of the numerical checks

The derivative-route and lag-correlation-support programs evaluate the paper's
QPAC inequalities in ordinary double precision. They are numerical
support-class checks, not directed-rounding interval certificates. 

The dual-polynomial experiment is a computational illustration of the
finite-harmonic lifted surrogate. 


## Reproducibility notes

- All physical parameters and source coefficients are defined inside the
  corresponding MATLAB functions.
- No generated file is required as an input to another experiment.
- Numerical values in the paper are rounded for presentation.
- Compare outputs using the reported tolerances rather than bitwise equality.
- Solver iteration counts and the last reported digits can depend on the
  MATLAB, BLAS, CVX, and SDPT3 versions.

## Citing the code
 Cite both the paper and the archived
software release.


## Rights and permitted use

Copyright (c) 2026 Sajad Daei. All rights reserved.

This code was written by Sajad Daei and is provided as supplementary
reproducibility material for

> S. Daei, G. Fodor, and M. Skoglund,  
> *A Mathematical Theory of Near-Field Super-Resolution*.

Permission is granted to download, inspect, and run the code for academic,
non-commercial peer review, verification, and reproducibility purposes.

Any publication, preprint, report, presentation, or software work that uses or
adapts this code, or uses numerical output derived from it, must cite the paper
and the archived software release.

Redistribution, modification, incorporation into other software, commercial
use, or public release of modified versions is not permitted without prior
written permission from Sajad Daei.

This code is provided “as is,” without warranty of any kind, express or
implied.
