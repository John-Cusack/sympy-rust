"""Equivalence corpus: Rust kernel vs sympy (exact rational equality).

Every case asserts bit-exact `QQ` agreement with upstream
(`densearith`/`Poly`), not just "no crash". Seeded RNG keeps the corpus
deterministic; edge rows (zero, constants, cancellation) are explicit.
"""

import random

import pytest
from sympy import QQ, Poly, expand
from sympy.abc import x, y, z
from sympy.polys.densearith import dup_add, dup_mul, dup_sub
from sympy.polys.polytools import gcd as poly_gcd_sympy

import sympy_poly
from sympy_poly import dmp_mul, dup_add as r_add
from sympy_poly import dup_gcd as r_gcd
from sympy_poly import dup_mul as r_mul
from sympy_poly import dup_sub as r_sub
from sympy_poly import poly_gcd, poly_mul

RNG = random.Random(20261017)


def rand_dup(deg, small=True):
    coeffs = []
    for _ in range(deg + 1):
        n = RNG.randint(-9, 9) if small else RNG.randint(-999, 999)
        d = RNG.randint(1, 9) if small else RNG.randint(1, 999)
        coeffs.append(QQ(n, d))
    if all(c == 0 for c in coeffs):
        coeffs[0] = QQ(1, 2)
    return coeffs


def test_dup_mul_matches_densearith():
    for deg in (0, 1, 2, 5, 12, 30):
        f, g = rand_dup(deg), rand_dup(deg)
        assert r_mul(f, g) == dup_mul(f, g, QQ)


def test_dup_mul_matches_poly():
    f = Poly(sum(c * x**i for i, c in enumerate(reversed(rand_dup(7)))), x, domain=QQ)
    g = Poly(sum(c * x**i for i, c in enumerate(reversed(rand_dup(7)))), x, domain=QQ)
    assert poly_mul(f, g) == f * g


def test_dup_add_sub_match():
    for deg in (0, 3, 15):
        f, g = rand_dup(deg), rand_dup(deg)
        assert r_add(f, g) == dup_add(f, g, QQ)
        assert r_sub(f, g) == dup_sub(f, g, QQ)


def test_dup_gcd_monic_associate():
    cases = [
        ([QQ(-1), QQ(0), QQ(1)], [QQ(-1), QQ(1)]),  # x^2-1, x-1
        ([QQ(2), QQ(2)], [QQ(3), QQ(3)]),  # associates -> x+1
        ([QQ(1), QQ(0), QQ(1)], [QQ(1), QQ(1)]),  # coprime -> 1
        ([], [QQ(2), QQ(2)]),  # gcd(0, f)
        ([], []),  # gcd(0, 0)
    ]
    for f, g in cases:
        assert r_gcd(f, g) == _monic_gcd_sympy(f, g)
    for _ in range(10):
        a = rand_dup(RNG.randint(1, 6))
        b = rand_dup(RNG.randint(1, 6))
        assert r_gcd(a, b) == _monic_gcd_sympy(a, b)

def _monic_gcd_sympy(f, g):
    from sympy.polys.densetools import dup_monic
    from sympy.polys.euclidtools import dup_gcd

    return dup_monic(dup_gcd(f, g, QQ), QQ)


def test_dup_gcd_matches_poly_gcd():
    a = Poly((x**2 - 1) * (x + 2), x, domain=QQ)
    b = Poly((x - 1) * (x + 2), x, domain=QQ)
    assert poly_gcd(a, b) == poly_gcd_sympy(a, b)


def test_dmp_mul_bivariate_matches_poly():
    e1 = expand((x + 2 * y + 1) ** 4)
    e2 = expand((x - y + 3) ** 3)
    p1, p2 = Poly(e1, x, y, domain=QQ), Poly(e2, x, y, domain=QQ)
    assert poly_mul(p1, p2) == p1 * p2


def test_dmp_mul_trivariate_matches_poly():
    e1 = expand((x + y + z + 1) ** 3)
    e2 = expand((x - y + 2 * z) ** 2)
    p1, p2 = Poly(e1, x, y, z, domain=QQ), Poly(e2, x, y, z, domain=QQ)
    assert poly_mul(p1, p2) == p1 * p2


def test_zero_and_constant_identities():
    one, zero = [QQ(1)], []
    f = rand_dup(5)
    assert r_mul(f, one) == f
    assert r_mul(f, zero) == []
    assert r_add(f, zero) == f
    assert r_sub(f, f) == []
    assert poly_mul(Poly(0, x, domain=QQ), Poly(x + 1, x, domain=QQ)) == 0
def test_zero_denominator_rejected():
    from sympy_poly import _core

    with pytest.raises(ZeroDivisionError):
        _core.dup_mul([(1, 1)], [(1, 0)])
