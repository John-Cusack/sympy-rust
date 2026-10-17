"""`sympy_poly`: Python face of the polys-domain kernel port.

Thin converters over the Rust `_core` extension. Coefficient lists use
sympy's dense order (high-to-low per level); rationals cross the seam as
`(num, den)` int pairs and come back as `(num_str, den_str)`.
"""

from sympy import QQ

from sympy_poly._core import dmp_mul as _dmp_mul
from sympy_poly._core import dup_add as _dup_add
from sympy_poly._core import dup_gcd as _dup_gcd
from sympy_poly._core import dup_mul as _dup_mul
from sympy_poly._core import dup_sub as _dup_sub

__all__ = ["dup_mul", "dup_add", "dup_sub", "dup_gcd", "dmp_mul", "poly_mul", "poly_gcd"]


def _to_pairs(coeffs):
    """`QQ` (or int) coeffs (high-to-low) -> `(num, den)` int pairs."""
    return [(int(c.numerator), int(c.denominator)) for c in coeffs]


def _from_pairs(pairs):
    """`(num_str, den_str)` pairs -> `QQ` coeffs (high-to-low)."""
    return [QQ(int(n), int(d)) for n, d in pairs]


def dup_mul(f, g):
    """Dense univariate multiply over QQ; lists of QQ, high-to-low."""
    return _from_pairs(_dup_mul(_to_pairs(f), _to_pairs(g)))


def dup_add(f, g):
    """Dense univariate add over QQ; lists of QQ, high-to-low."""
    return _from_pairs(_dup_add(_to_pairs(f), _to_pairs(g)))


def dup_sub(f, g):
    """Dense univariate subtract over QQ; lists of QQ, high-to-low."""
    return _from_pairs(_dup_sub(_to_pairs(f), _to_pairs(g)))


def dup_gcd(f, g):
    """Dense univariate GCD over QQ, monic; lists of QQ, high-to-low."""
    return _from_pairs(_dup_gcd(_to_pairs(f), _to_pairs(g)))


def _nested_to_pairs(nested):
    if isinstance(nested, list) and nested and isinstance(nested[0], list):
        return [_nested_to_pairs(c) for c in nested]
    return _to_pairs(nested)


def _nested_from_pairs(nested):
    if isinstance(nested, list) and nested and isinstance(nested[0], list):
        return [_nested_from_pairs(c) for c in nested]
    return _from_pairs(nested)


def dmp_mul(a, b, levels):
    """Recursive dense multiply over QQ; nested lists, high-to-low."""
    return _nested_from_pairs(_dmp_mul(_nested_to_pairs(a), _nested_to_pairs(b), levels))


def poly_mul(a, b):
    """`sympy.Poly` multiply over QQ via the Rust kernel."""
    from sympy import Poly

    if not isinstance(a, Poly):
        a = Poly(a)
    if not isinstance(b, Poly):
        b = Poly(b)
    gens, domain = a.gens, a.domain
    assert b.gens == gens and a.domain == domain == QQ, "QQ polys over one generator set only"
    if len(gens) == 1:
        coeffs = dup_mul(a.all_coeffs(), b.all_coeffs())
        return Poly(coeffs, *gens, domain=QQ)
    nested_a = _poly_to_nested(a, len(gens))
    nested_b = _poly_to_nested(b, len(gens))
    return _nested_to_poly(dmp_mul(nested_a, nested_b, len(gens)), gens)

def _poly_to_nested(p, levels):
    """`Poly` -> dense nested lists (high-to-low per level) via `as_dict`."""
    from sympy import QQ

    d = p.as_dict()
    widths = [0] * levels
    for exps in d:
        for i, e in enumerate(exps):
            widths[i] = max(widths[i], e)

    def build(level, prefix):
        w = widths[level]
        if level == levels - 1:
            return [d.get(prefix + (e,), QQ(0)) for e in range(w, -1, -1)]
        return [build(level + 1, prefix + (e,)) for e in range(w, -1, -1)]

    return build(0, ())


def _nested_to_poly(nested, gens):
    """Dense nested lists (any shape) -> `Poly` via an exponent dict."""
    from sympy import Poly, QQ

    d = {}

    def walk(node, level, prefix):
        last = len(gens) - 1
        w = len(node) - 1
        for i, c in enumerate(node):
            exp = prefix + (w - i,)
            if level == last:
                if c != 0:
                    d[exp] = c
            else:
                walk(c, level + 1, exp)

    walk(nested, 0, ())
    if not d:
        return Poly(QQ(0), *gens, domain=QQ)
    return Poly(d, *gens, domain=QQ)


def poly_gcd(a, b):
    """Univariate `sympy.Poly` GCD over QQ via the Rust kernel (monic)."""
    from sympy import Poly

    if not isinstance(a, Poly):
        a = Poly(a)
    if not isinstance(b, Poly):
        b = Poly(b)
    assert len(a.gens) == 1 and len(b.gens) == 1, "univariate only"
    return Poly(dup_gcd(a.all_coeffs(), b.all_coeffs()), a.gens[0], domain=QQ)
