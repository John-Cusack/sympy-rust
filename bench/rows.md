# bench (measured-only)

machine: x86_64, python 3.12.3, sympy 1.14.0
wall = best-of-3 of the kernel op; maxrss = process peak (setup included, same shape both sides).

| workload | sympy wall (s) | rust wall (s) | speedup | sympy RSS (MB) | rust RSS (MB) |
|---|---|---|---|---|---|
| dup_mul_800 | 0.212 | 0.042 | 5.06x | 46.3 | 46.9 |
| dmp_mul_3var | 0.112 | 0.013 | 8.75x | 48.6 | 48.9 |
| dup_gcd_70 | 0.001 | 0.008 | 0.17x | 49.0 | 48.9 |
