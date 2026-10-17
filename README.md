# sympy-rust

Rust port of the [`sympy/sympy`](https://github.com/sympy/sympy)
polys-domain kernel: dense polynomial arithmetic over `QQ` plus the `QQ`
domain operations. See `NOTICE` for upstream attribution.

## Scope

Kernel-scoped, not a mirror bid. SymPy is 497K LOC of pure Python; the
rustsmith headroom probe graded the two candidate kernels:

| kernel | score | verdict |
|---|---|---|
| polys-domains (`dup_mul` + `pythonmpq` QQ ops) | 0.85 | rewrite |
| simplify-core (9-module dispatch) | 0.49 | upstream_only — out of scope |

Port seam: `dup_mul`/`dup_add`/`dup_sub`/`dup_gcd` (univariate dense over
QQ) + recursive `dmp_mul` (multivariate) + `QQ` exact rationals. GCD
output is monic (unit-normalized, as upstream). Expand only with measured
wins.

## Layout

- `sympy-poly-core/` — reusable Rust core, no Python dependency.
- `src/lib.rs` — thin PyO3 binding (`sympy_poly._core`).
- `sympy_poly/` — Python converters (`Poly` ↔ dense lists).
- `tests/test_equivalence.py` — exact-equality corpus vs sympy.
- `bench/` — two-track wall + peak-RSS benches, measured-only.

## Install (from source)

```sh
pip install maturin
maturin develop   # inside a virtualenv
pytest tests/
python bench/bench_poly.py
```

## Status

Kernel slice v0.1.0, measured on x86_64 / python 3.12.3 / sympy 1.14.0
(`bench/results.json`; wall = best-of-3, maxrss = process peak):

| workload | sympy wall (s) | rust wall (s) | speedup | sympy RSS (MB) | rust RSS (MB) |
|---|---|---|---|---|---|
| dup_mul_800 (univariate deg-800 dense QQ) | 0.212 | 0.042 | 5.06x | 46.3 | 46.9 |
| dmp_mul_3var (`(x+y+z+1)^12` squared) | 0.112 | 0.013 | 8.75x | 48.6 | 48.9 |
| dup_gcd_70 (heuristic-friendly shape) | 0.001 | 0.008 | 0.17x, slower | 49.0 | 48.9 |

Accuracy: `tests/test_equivalence.py` 9/9 exact vs sympy; core unit
tests 14/14 (`cargo test -p sympy-poly-core`); upstream kernel subset
(`test_densearith` + `test_densetools` + `test_euclidtools`) 105 passed
on origin as baseline.

Ranked follow-ups: heuristic (evaluation/interpolation) GCD for
heuristic-friendly shapes — the PRS here is correct but ties-or-loses
where sympy's heuristic fires; multivariate GCD; `factor` stays upstream.
Publish (TestPyPI/PyPI/crates.io) is human-gated and pending.
