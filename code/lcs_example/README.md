# Lag-correlation support (LCS) QPAC numerical test

Run from this folder:

```matlab
result = run_lcs_example_qpac();
```

The script implements the same-bearing, different-range LCS example. The support-support LCS calculation uses analytic magnitude/sign bounds for the structured `K`, `H`, `dK`, and `dH` coefficient products. The near-region calculation uses the full Taylor segment tube and physical first angular derivatives in `C1`.

