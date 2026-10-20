# Release spec — `sympy-poly-rust` 0.1.0

Release lane spec for the kernel-scoped SymPy polys port. It records what
ships, how it is verified, and what it does not claim. The lane-equivalent
config is `release.toml`; attribution is `NOTICE`; this port's own license
is `LICENSE` (MIT).

## Scope

- Python distribution `sympy-poly-rust` (import name `sympy_poly`), version
  `0.1.0`.
- Reusable Rust core crate `sympy-poly-core` 0.1.0 on crates.io.
- Kernel seam only: dense univariate/multivariate arithmetic over `QQ`
  (`dup_mul`, `dup_add`, `dup_sub`, `dup_gcd`, `dmp_mul`) plus the `QQ`
  domain operations. `simplify` and `factor` stay upstream.

This is **not** a drop-in replacement for `sympy`: the upstream package
name and import stay with upstream, and the port covers one kernel seam.
The published distribution is named accordingly.

## Artifacts

| Artifact | Content |
|---|---|
| wheel | `sympy_poly` package with the compiled `sympy_poly._core` extension (PyO3) |
| sdist | source tree incl. `sympy-poly-core/` (path dependency), `NOTICE`, `LICENSE`, `release.toml` |
| crate | `sympy-poly-core` 0.1.0 (rlib, no Python dependency) |

## Verification battery (all green before any upload)

1. `cargo test -p sympy-poly-core` — the core unit tests Rust consumers use.
2. A throwaway downstream crate depending on the core by path, built with
   no Python anywhere — the core is reusable as a plain `rlib`.
3. `cargo package` on the core — what crates.io receives.
4. `maturin build --release` — a release-profile wheel (never a develop
   build).
5. `maturin sdist` — the source distribution, carrying the core.
6. Clean venv, wheel installed from disk: dist name/version match, and every
   `smoke_exprs` entry from `release.toml` holds.
7. Sdist unpacked, wheel rebuilt from the unpacked sdist (proves the core
   ships inside), and that rebuilt wheel installed into a second clean venv
   with the same smoke vectors — the sdist is verified, not assumed.
8. The port's own parity corpus (`tests/test_equivalence.py`, 9/9 exact vs
   sympy) re-run against the installed wheel.

## Registry lanes

1. TestPyPI first; install from the TestPyPI index into a clean venv and
   re-run the smoke vectors (the index, not the local file, is proven).
2. Production PyPI: the same wheel and sdist, same SHA-256, promoted only
   after the TestPyPI success. The publisher account carries a PyPI
   new-project throttle; a `429` on the first production create is parked
   (no retry), with the exact error recorded.
3. crates.io: `sympy-poly-core` 0.1.0, one attempt.

## What this release does not claim

- No speedup outside the measured workloads in the README table
  (`bench/results.json`, measured 2026-10-17).
- `dup_gcd` is slower than sympy on heuristic-friendly shapes (0.17x,
  measured); the PRS is correct but ties-or-loses where sympy's heuristic
  fires. No claim is made there.
- No claim that the port replaces `sympy` outside the kernel seam.
