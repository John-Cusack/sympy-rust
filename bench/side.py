"""One benchmark side: build inputs in setup, time the kernel op, print JSON.

Usage: side.py <sympy|rust> <dup_mul_800|dmp_mul_3var|dup_gcd_40>

`wall_secs` is best-of-3 of the op only. `maxrss_kb` is the process
lifetime peak (setup included, same shape both sides — compare
relatively, not absolutely).
"""

import json
import resource
import sys
import time

SIDE, WORKLOAD = sys.argv[1], sys.argv[2]


def build_dup_mul_800():
    import random
    from sympy import QQ

    rng = random.Random(77)
    f = [QQ(rng.randint(-9, 9), rng.randint(1, 9)) for _ in range(801)]
    g = [QQ(rng.randint(-9, 9), rng.randint(1, 9)) for _ in range(801)]
    return f, g


def build_dmp_mul_3var():
    from sympy import Poly, QQ, expand
    from sympy.abc import x, y, z

    p = Poly(expand((x + y + z + 1) ** 12), x, y, z, domain=QQ)
    return p, p


def build_dup_gcd_70():
    import random
    from sympy import QQ, Poly
    from sympy.abc import x

    rng = random.Random(99)

    def rand_poly(deg):
        return Poly(
            sum(QQ(rng.randint(-5, 5), rng.randint(1, 5)) * x**i for i in range(deg + 1)), x
        )

    # Heuristic-friendly shape: large shared factor; sympy fires its
    # evaluation/interpolation GCD here, our PRS does not (ranked follow-up).
    common = rand_poly(30)
    return common * rand_poly(60), common * rand_poly(70)


BUILDERS = {
    "dup_mul_800": build_dup_mul_800,
    "dmp_mul_3var": build_dmp_mul_3var,
    "dup_gcd_70": build_dup_gcd_70,
}


def op_sympy(workload, inputs):
    if workload == "dup_mul_800":
        from sympy import QQ
        from sympy.polys.densearith import dup_mul

        f, g = inputs
        return lambda: dup_mul(f, g, QQ)
    if workload == "dmp_mul_3var":
        p, _ = inputs
        return lambda: p * p
    from sympy.polys.polytools import gcd as pgcd

    a, b = inputs
    return lambda: pgcd(a, b)


def op_rust(workload, inputs):
    import sympy_poly

    if workload == "dup_mul_800":
        f, g = inputs
        return lambda: sympy_poly.dup_mul(f, g)
    if workload == "dmp_mul_3var":
        from sympy_poly import _poly_to_nested

        p, _ = inputs
        gens = p.gens
        nested = _poly_to_nested(p, len(gens))
        return lambda: sympy_poly.dmp_mul(nested, nested, len(gens))
    a, b = inputs
    return lambda: sympy_poly.poly_gcd(a, b)


def main():
    inputs = BUILDERS[WORKLOAD]()
    run = op_rust(WORKLOAD, inputs) if SIDE == "rust" else op_sympy(WORKLOAD, inputs)
    run()  # warmup (caches, allocator)
    best = min(_time(run) for _ in range(3))
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps({"side": SIDE, "workload": WORKLOAD, "wall_secs": best, "maxrss_kb": peak_kb}))


def _time(fn):
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


main()
