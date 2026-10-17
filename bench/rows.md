# bench (measured-only)

machine: x86_64, python 3.12.3, sympy 1.14.0
wall = best-of-3 of the kernel op; maxrss = process peak (setup included, same shape both sides).

| workload | sympy wall (s) | rust wall (s) | speedup | sympy RSS (MB) | rust RSS (MB) |
|---|---|---|---|---|---|
| dup_mul_800 | 0.213 | 0.040 | 5.31x | 46.5 | 47.2 |
| dmp_mul_3var | 0.108 | 0.012 | 8.85x | 48.6 | 49.9 |
| dup_gcd_70 | 0.002 | 0.008 | 0.20x | 46.7 | 48.4 |
