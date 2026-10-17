//! Exact integer polynomial arithmetic backing the QQ GCD.
//!
//! Mirrors sympy's `dup_prem` (`densearith`) and the subresultant PRS
//! (`dup_inner_subresultants`, `dup_rr_prs_gcd` in `euclidtools`): the QQ
//! GCD clears denominators, runs the fraction-free PRS over ZZ, and
//! returns the monic associate. Naive field-Euclid blows coefficients up
//! exponentially; the PRS keeps them polynomial-size.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

/// Dense univariate over ZZ, low-degree first; empty = zero.
pub type ZPoly = Vec<BigInt>;

pub fn strip(mut f: ZPoly) -> ZPoly {
    while f.last().is_some_and(|c| c.is_zero()) {
        f.pop();
    }
    f
}

/// Degree, or `-1` for the zero polynomial.
pub fn degree(f: &[BigInt]) -> i64 {
    f.len() as i64 - 1
}

fn lc(f: &[BigInt]) -> BigInt {
    f.last().cloned().unwrap_or_else(BigInt::zero)
}

fn mul_ground(f: &[BigInt], c: &BigInt) -> ZPoly {
    if c.is_zero() {
        return Vec::new();
    }
    strip(f.iter().map(|a| a * c).collect())
}

/// Exact ground division (exact by construction at every call site).
fn quo_ground(f: &[BigInt], c: &BigInt) -> ZPoly {
    assert!(!c.is_zero(), "ground division by zero");
    strip(
        f.iter()
            .map(|a| {
                let (q, r) = a.div_rem(c);
                assert!(r.is_zero(), "inexact ground division");
                q
            })
            .collect(),
    )
}

/// Pseudo-remainder of `f` by `g` (mirrors `dup_prem`).
pub fn prem(f: &[BigInt], g: &[BigInt]) -> ZPoly {
    assert!(!g.is_empty(), "pseudo-division by zero");
    let (df, dg) = (degree(f), degree(g));
    if df < dg {
        return f.to_vec();
    }
    let mut r = f.to_vec();
    let mut dr = df;
    let mut n = df - dg + 1;
    let lc_g = lc(g);
    loop {
        let lc_r = lc(&r);
        let j = (dr - dg) as usize;
        n -= 1;
        let mut next = vec![BigInt::zero(); r.len().max(j + g.len())];
        for (i, c) in r.iter().enumerate() {
            next[i] += c * &lc_g;
        }
        for (i, c) in g.iter().enumerate() {
            next[i + j] -= c * &lc_r;
        }
        r = strip(next);
        let prev = dr;
        dr = degree(&r);
        if dr < dg {
            break;
        }
        assert!(dr < prev, "prem failed to reduce");
    }
    mul_ground(&r, &lc_g.pow(n.max(0) as u32))
}

/// Content (gcd of coefficients, non-negative) and primitive part.
pub fn primitive(f: &[BigInt]) -> (BigInt, ZPoly) {
    let mut c = BigInt::zero();
    for a in f {
        c = c.gcd(&a.abs());
        if c.is_one() {
            break;
        }
    }
    if c.is_zero() {
        return (BigInt::zero(), Vec::new());
    }
    (c.clone(), quo_ground(f, &c))
}

/// Last (constant-or-higher) element of the subresultant PRS, mirroring
/// `dup_inner_subresultants` (low-to-high adaptation; same arithmetic).
pub fn subresultant_last(f0: &[BigInt], g0: &[BigInt]) -> ZPoly {
    let (mut f, mut g) = if degree(f0) >= degree(g0) {
        (f0.to_vec(), g0.to_vec())
    } else {
        (g0.to_vec(), f0.to_vec())
    };
    if f.is_empty() {
        return Vec::new();
    }
    if g.is_empty() {
        return f;
    }
    let (n, mut m) = (degree(&f), degree(&g));
    let mut d = n - m;
    let mut b = BigInt::from(-1).pow((d + 1).max(0) as u32);
    let mut h = mul_ground(&prem(&f, &g), &b);
    let mut lc_cur = lc(&g);
    let mut c = -lc_cur.pow(d.max(0) as u32);
    while !h.is_empty() {
        let k = degree(&h);
        f = g;
        g = h;
        d = m - k;
        m = k;
        b = -&lc_cur * c.pow(d.max(0) as u32);
        h = quo_ground(&prem(&f, &g), &b);
        lc_cur = lc(&g);
        if d > 1 {
            let q = c.pow((d - 1) as u32);
            c = (-&lc_cur).pow(d as u32).div_rem(&q).0;
        } else {
            c = -lc_cur.clone();
        }
    }
    g
}

/// ZZ GCD: content gcd times the primitive subresultant, positive lc.
pub fn gcd(f: &[BigInt], g: &[BigInt]) -> ZPoly {
    if f.is_empty() && g.is_empty() {
        return Vec::new();
    }
    if f.is_empty() || g.is_empty() {
        let h = if f.is_empty() { g } else { f };
        let (_, mut p) = primitive(h);
        if lc(&p).is_negative() {
            p = p.into_iter().map(|a| -a).collect();
        }
        return p;
    }
    let (cf, pf) = primitive(f);
    let (cg, pg) = primitive(g);
    let c = cf.gcd(&cg);
    let (_, mut h) = primitive(&subresultant_last(&pf, &pg));
    if lc(&h).is_negative() {
        h = h.into_iter().map(|a| -a).collect();
    }
    mul_ground(&h, &c)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn z(ns: &[i64]) -> ZPoly {
        ns.iter().map(|&n| BigInt::from(n)).collect()
    }

    #[test]
    fn prem_matches_upstream() {
        // dup_prem(x^2 + 1, 2x - 4) = 20
        assert_eq!(prem(&z(&[1, 0, 1]), &z(&[-4, 2])), z(&[20]));
        // deg f < deg g returns f
        assert_eq!(prem(&z(&[1, 1]), &z(&[1, 0, 1])), z(&[1, 1]));
    }

    #[test]
    fn subresultant_last_matches_upstream() {
        // dup_subresultants(x^2 + 1, x^2 - 1) ends in -2
        assert_eq!(
            subresultant_last(&z(&[1, 0, 1]), &z(&[-1, 0, 1])),
            z(&[-2])
        );
    }

    #[test]
    fn gcd_is_primitive_positive() {
        // gcd(x^2 - 1, x^2 - 3x + 2) = x - 1
        assert_eq!(
            gcd(&z(&[-1, 0, 1]), &z(&[2, -3, 1])),
            z(&[-1, 1])
        );
        // scaled inputs: gcd(2x + 2, 3x + 3) = x + 1
        assert_eq!(gcd(&z(&[2, 2]), &z(&[3, 3])), z(&[1, 1]));
        // coprime
        assert_eq!(gcd(&z(&[1, 0, 1]), &z(&[1, 1])), z(&[1]));
        assert!(gcd(&[], &[]).is_empty());
        assert_eq!(gcd(&[], &z(&[2, 2])), z(&[1, 1]));
    }
}
