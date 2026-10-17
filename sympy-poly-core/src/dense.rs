//! Dense polynomial arithmetic over `QQ`, mirroring sympy's
//! `polys/densearith.py` (`dup_*`) and the recursive `dmp_*` layer.
//!
//! Conventions (same as upstream dense polys):
//!
//! * univariate coefficient vecs are stored **low-degree first**
//!   (`v[i]` is the `x^i` coefficient); the binding converts to/from
//!   sympy's high-to-low lists at the seam.
//! * the zero polynomial is the empty vec; [`strip`] restores that form.
//! * `dmp` recursion nests [`Dmp::Poly`] levels outermost-first with the
//!   same low-to-high order per level; [`Dmp::Coeff`] is a depth-0 leaf.

use crate::qq::Qq;

/// Recursive dense multivariate polynomial over `QQ`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Dmp {
    /// Depth-0 leaf: a bare domain element.
    Coeff(Qq),
    /// One dense level, low-degree first; empty = zero at this depth.
    Poly(Vec<Dmp>),
}

impl Dmp {
    pub fn zero(depth: u32) -> Self {
        if depth == 0 {
            Self::Coeff(Qq::zero())
        } else {
            Self::Poly(Vec::new())
        }
    }

    pub fn is_zero(&self) -> bool {
        match self {
            Self::Coeff(q) => q.is_zero(),
            Self::Poly(cs) => cs.iter().all(Dmp::is_zero),
        }
    }

    /// Drop leading zero children, restoring canonical form.
    pub fn strip(self) -> Self {
        match self {
            Self::Coeff(q) => Self::Coeff(q),
            Self::Poly(mut cs) => {
                while cs.last().is_some_and(Dmp::is_zero) {
                    cs.pop();
                }
                Self::Poly(cs)
            }
        }
    }

    /// Degree in the outermost variable (`-1` convention avoided: zero
    /// maps to `None`, matching "no terms").
    pub fn degree(&self) -> Option<usize> {
        match self {
            Self::Coeff(q) => {
                if q.is_zero() {
                    None
                } else {
                    Some(0)
                }
            }
            Self::Poly(cs) => {
                if cs.is_empty() {
                    None
                } else {
                    Some(cs.len() - 1)
                }
            }
        }
    }
}

/// Univariate dense multiply (schoolbook, zero-skipping like `dup_mul`).
pub fn dup_mul(f: &[Qq], g: &[Qq]) -> Vec<Qq> {
    if f.is_empty() || g.is_empty() {
        return Vec::new();
    }
    let mut h = vec![Qq::zero(); f.len() + g.len() - 1];
    for (i, a) in f.iter().enumerate() {
        if a.is_zero() {
            continue;
        }
        for (j, b) in g.iter().enumerate() {
            if b.is_zero() {
                continue;
            }
            h[i + j] = &h[i + j] + &(a * b);
        }
    }
    strip_vec(h)
}

/// Univariate dense add.
pub fn dup_add(f: &[Qq], g: &[Qq]) -> Vec<Qq> {
    let n = f.len().max(g.len());
    let zero = Qq::zero();
    let mut h = Vec::with_capacity(n);
    for i in 0..n {
        let a = f.get(i).unwrap_or(&zero);
        let b = g.get(i).unwrap_or(&zero);
        h.push(a + b);
    }
    strip_vec(h)
}

/// Univariate dense subtract.
pub fn dup_sub(f: &[Qq], g: &[Qq]) -> Vec<Qq> {
    let n = f.len().max(g.len());
    let zero = Qq::zero();
    let mut h = Vec::with_capacity(n);
    for i in 0..n {
        let a = f.get(i).unwrap_or(&zero);
        let b = g.get(i).unwrap_or(&zero);
        h.push(a - b);
    }
    strip_vec(h)
}

/// Univariate dense negate.
pub fn dup_neg(f: &[Qq]) -> Vec<Qq> {
    strip_vec(f.iter().map(|a| -a).collect())
}

fn strip_vec(mut f: Vec<Qq>) -> Vec<Qq> {
    while f.last().is_some_and(Qq::is_zero) {
        f.pop();
    }
    f
}

/// Exact univariate division over a field: `a = q*b + r`, `deg r < deg b`.
/// Returns `None` on division by zero.
pub fn dup_div(a: &[Qq], b: &[Qq]) -> Option<(Vec<Qq>, Vec<Qq>)> {
    let b = strip_vec(b.to_vec());
    if b.is_empty() {
        return None;
    }
    let mut r = strip_vec(a.to_vec());
    let db = b.len() - 1;
    let lc = b.last().expect("non-empty divisor");
    let mut q = vec![Qq::zero(); r.len().saturating_sub(b.len()) + 1];
    while r.len() >= b.len() && !r.is_empty() {
        let dr = r.len() - 1;
        let coeff = r[dr].exact_div(lc)?;
        let shift = dr - db;
        q[shift] = &q[shift] + &coeff;
        for (j, bj) in b.iter().enumerate() {
            let k = shift + j;
            r[k] = &r[k] - &(&coeff * bj);
        }
        r = strip_vec(r);
    }
    Some((strip_vec(q), r))
}

/// Univariate GCD over `QQ`, monic output (unit-normalized like the
/// domain GCD: associates compare equal after making monic).
/// Zero convention matches upstream: `gcd(0, 0) = 0`, `gcd(0, f)` is the
/// monic associate of `f`.
///
/// Route: clear denominators, run the fraction-free subresultant PRS
/// over ZZ (mirrors `dup_rr_prs_gcd`), return the monic associate.
/// Naive field-Euclid blows coefficients up exponentially; the PRS
/// keeps them polynomial-size.
pub fn dup_gcd(a: &[Qq], b: &[Qq]) -> Vec<Qq> {
    let a = strip_vec(a.to_vec());
    let b = strip_vec(b.to_vec());
    if a.is_empty() {
        return dup_monic(&b);
    }
    if b.is_empty() {
        return dup_monic(&a);
    }
    if a.len() == 1 || b.len() == 1 {
        // Nonzero constants are units over QQ: the gcd is 1.
        return vec![Qq::one()];
    }
    let h = crate::zz::gcd(&to_zz(&a), &to_zz(&b));
    let q: Vec<Qq> = h
        .iter()
        .map(|c| Qq::new(c.clone(), num_bigint::BigInt::from(1)).expect("unit denominator"))
        .collect();
    dup_monic(&q)
}

/// Clear denominators: scale to a ZZ polynomial (low-to-high).
fn to_zz(f: &[Qq]) -> crate::zz::ZPoly {
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::One;
    let mut l = BigInt::one();
    for c in f {
        l = l.lcm(&c.denom());
    }
    f.iter()
        .map(|c| c.numer() * (l.clone() / c.denom()))
        .collect()
}

/// Scale to monic (leading coefficient 1); zero stays zero.
pub fn dup_monic(f: &[Qq]) -> Vec<Qq> {
    let f = strip_vec(f.to_vec());
    match f.last() {
        None => Vec::new(),
        Some(lc) => {
            if lc.is_one() {
                return f;
            }
            f.iter()
                .map(|c| c.exact_div(lc).expect("field division is exact"))
                .collect()
        }
    }
}

/// Recursive dense multiply at `depth` remaining levels (mirrors `dmp_mul`).
pub fn dmp_mul(a: &Dmp, b: &Dmp, depth: u32) -> Dmp {
    if a.is_zero() || b.is_zero() {
        return Dmp::zero(depth);
    }
    if depth == 0 {
        match (a, b) {
            (Dmp::Coeff(x), Dmp::Coeff(y)) => return Dmp::Coeff(x * y),
            _ => panic!("depth-0 dmp_mul needs Coeff leaves"),
        }
    }
    let (Dmp::Poly(f), Dmp::Poly(g)) = (a, b) else {
        panic!("depth-{depth} dmp_mul needs Poly levels");
    };
    let mut h: Vec<Dmp> = (0..f.len() + g.len() - 1)
        .map(|_| Dmp::zero(depth - 1))
        .collect();
    for (i, x) in f.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in g.iter().enumerate() {
            if y.is_zero() {
                continue;
            }
            let p = dmp_mul(x, y, depth - 1);
            h[i + j] = dmp_add(&h[i + j], &p, depth - 1);
        }
    }
    Dmp::Poly(h).strip()
}

/// Recursive dense add at `depth` remaining levels.
pub fn dmp_add(a: &Dmp, b: &Dmp, depth: u32) -> Dmp {
    if depth == 0 {
        match (a, b) {
            (Dmp::Coeff(x), Dmp::Coeff(y)) => return Dmp::Coeff(x + y),
            _ => panic!("depth-0 dmp_add needs Coeff leaves"),
        }
    }
    let (Dmp::Poly(f), Dmp::Poly(g)) = (a, b) else {
        panic!("depth-{depth} dmp_add needs Poly levels");
    };
    let n = f.len().max(g.len());
    let z = Dmp::zero(depth - 1);
    let mut h = Vec::with_capacity(n);
    for i in 0..n {
        let x = f.get(i).unwrap_or(&z);
        let y = g.get(i).unwrap_or(&z);
        h.push(dmp_add(x, y, depth - 1));
    }
    Dmp::Poly(h).strip()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    fn qq(n: i64, d: i64) -> Qq {
        Qq::new(BigInt::from(n), BigInt::from(d)).unwrap()
    }

    /// low-to-high helper: ints as QQ.
    fn v(ns: &[i64]) -> Vec<Qq> {
        ns.iter().map(|&n| Qq::from_int(n)).collect()
    }

    #[test]
    fn mul_matches_schoolbook() {
        // (1 + 2x)(3 + 4x) = 3 + 10x + 8x^2
        assert_eq!(dup_mul(&v(&[1, 2]), &v(&[3, 4])), v(&[3, 10, 8]));
        // rationals: (1/2 + x)(2x) = x + 2x^2
        assert_eq!(
            dup_mul(&[qq(1, 2), Qq::one()], &[Qq::zero(), qq(2, 1)]),
            vec![Qq::zero(), Qq::one(), qq(2, 1)]
        );
    }

    #[test]
    fn zero_and_empty_identities() {
        assert!(dup_mul(&v(&[1, 2]), &[]).is_empty());
        assert!(dup_mul(&[], &v(&[1])).is_empty());
        // zero polynomial annihilates even with trailing-zero junk
        assert!(dup_mul(&[Qq::zero(), Qq::zero()], &v(&[1, 1])).is_empty());
    }

    #[test]
    fn add_sub_neg_strip() {
        assert_eq!(dup_add(&v(&[1, 2]), &v(&[3])), v(&[4, 2]));
        assert_eq!(dup_sub(&v(&[1, 2]), &v(&[1, 2])), Vec::new());
        assert_eq!(dup_neg(&v(&[1, -2])), v(&[-1, 2]));
        // degree mismatch + cancellation strips
        assert_eq!(dup_add(&v(&[1, 1]), &v(&[-1])), v(&[0, 1]));
    }

    #[test]
    fn div_is_exact_quotient_and_remainder() {
        // (x^2 - 1) / (x - 1) = (x + 1) rem 0
        let (q, r) = dup_div(&v(&[-1, 0, 1]), &v(&[-1, 1])).unwrap();
        assert_eq!(q, v(&[1, 1]));
        assert!(r.is_empty());
        // (x^2 + 1) / (x + 1): rem 2
        let (q, r) = dup_div(&v(&[1, 0, 1]), &v(&[1, 1])).unwrap();
        assert_eq!(q, v(&[-1, 1]));
        assert_eq!(r, v(&[2]));
        assert!(dup_div(&v(&[1]), &[]).is_none());
    }

    #[test]
    fn gcd_is_monic_associate() {
        // gcd(x^2 - 1, x - 1) = x - 1
        assert_eq!(dup_gcd(&v(&[-1, 0, 1]), &v(&[-1, 1])), v(&[-1, 1]));
        // scaled inputs still give the monic associate: gcd(2x+2, 3x+3) = x+1
        assert_eq!(dup_gcd(&v(&[2, 2]), &v(&[3, 3])), v(&[1, 1]));
        // coprime
        assert_eq!(dup_gcd(&v(&[1, 0, 1]), &v(&[1, 1])), v(&[1]));
        // gcd(0, f) = monic(f)
        assert_eq!(dup_gcd(&[], &v(&[2, 2])), v(&[1, 1]));
        assert!(dup_gcd(&[], &[]).is_empty());
    }

    #[test]
    fn dmp_mul_bivariate() {
        // (x + y)(x - y) = x^2 - y^2, outer var x, inner y, low-to-high.
        fn c(n: i64, d: i64) -> Dmp {
            Dmp::Coeff(qq(n, d))
        }
        fn p(cs: Vec<Dmp>) -> Dmp {
            Dmp::Poly(cs)
        }
        let xpy = p(vec![p(vec![c(0, 1), c(1, 1)]), p(vec![c(1, 1)])]);
        let xmy = p(vec![p(vec![c(0, 1), c(-1, 1)]), p(vec![c(1, 1)])]);
        let h = dmp_mul(&xpy, &xmy, 2);
        let expect = p(vec![
            p(vec![c(0, 1), c(0, 1), c(-1, 1)]), // -y^2
            p(vec![]),                            // 0*x (stripped to empty)
            p(vec![c(1, 1)]),                   // x^2
        ]);
        assert_eq!(h, expect.strip());
        // zero annihilates, one is identity
        assert!(dmp_mul(&xpy, &Dmp::zero(2), 2).is_zero());
    }
}
