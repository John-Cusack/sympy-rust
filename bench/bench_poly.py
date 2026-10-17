"""Two-track bench: sympy vs Rust kernel, wall time + peak RSS per side.

Each side runs in its own process (`bench/side.py`) so peak RSS is
attributable. Machine, versions, and raw figures land in
`bench/results.json`; the rendered table in `bench/rows.md`. All figures
measured on this box — no modeled speedups.
"""

import json
import platform
import subprocess
import sys
from pathlib import Path
HERE = Path(__file__).parent
WORKLOADS = ["dup_mul_800", "dmp_mul_3var", "dup_gcd_70"]
PYTHON = sys.executable


def run_side(side, workload):
    out = subprocess.run(
        [PYTHON, str(HERE / "side.py"), side, workload],
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(out.stdout)


def main():
    import sympy

    rows = []
    for wl in WORKLOADS:
        sym = run_side("sympy", wl)
        rst = run_side("rust", wl)
        rows.append(
            {
                "workload": wl,
                "sympy_wall_secs": sym["wall_secs"],
                "rust_wall_secs": rst["wall_secs"],
                "speedup": sym["wall_secs"] / rst["wall_secs"] if rst["wall_secs"] > 0 else None,
                "sympy_maxrss_kb": sym["maxrss_kb"],
                "rust_maxrss_kb": rst["maxrss_kb"],
            }
        )
    results = {
        "machine": platform.processor() or platform.machine(),
        "python": platform.python_version(),
        "sympy": sympy.__version__,
        "rows": rows,
    }
    (HERE / "results.json").write_text(json.dumps(results, indent=2) + "\n")

    lines = [
        "# bench (measured-only)",
        "",
        f"machine: {results['machine']}, python {results['python']}, sympy {results['sympy']}",
        "wall = best-of-3 of the kernel op; maxrss = process peak (setup included, same shape both sides).",
        "",
        "| workload | sympy wall (s) | rust wall (s) | speedup | sympy RSS (MB) | rust RSS (MB) |",
        "|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['workload']} | {r['sympy_wall_secs']:.3f} | {r['rust_wall_secs']:.3f} "
            f"| {r['speedup']:.2f}x | {r['sympy_maxrss_kb']/1024:.1f} | {r['rust_maxrss_kb']/1024:.1f} |"
        )
    (HERE / "rows.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


main()
