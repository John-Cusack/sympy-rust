//! `QQ`: exact rationals, mirroring sympy's `pythonmpq`.
//!
//! Representation (the lesson of the first bench: generic `BigInt` loses
//! 3x to C small-int fast paths on CAS-sized coefficients): normalized
//! `i64` pairs while values fit, `Ratio<BigInt>` past overflow. Every op
//! normalizes, like upstream.
//!
//! Invariant: numerator/denominator are coprime, denominator is strictly
//! positive, zero is `Small(0, 1)`. Division by zero is `None`, never a
//! panic: callers surface it as a Python `ZeroDivisionError`.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, ToPrimitive, Zero};

/// An exact rational (the `QQ` domain element).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Qq {
    Small(i64, i64),
    Big(Ratio<BigInt>),
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// Normalize an `i128` pair, downshifting to `Small` when it fits.
fn small(num: i128, den: i128) -> Qq {
    assert!(den != 0, "zero denominator");
    let (mut n, mut d) = (num, den);
    if d < 0 {
        n = -n;
        d = -d;
    }
    let g = gcd_u128(n.unsigned_abs(), d as u128) as i128;
    n /= g;
    d /= g;
    if i64::MIN as i128 <= n && n <= i64::MAX as i128 && d <= i64::MAX as i128 {
        Qq::Small(n as i64, d as i64)
    } else {
        Qq::Big(Ratio::new(BigInt::from(n), BigInt::from(d)))
    }
}

fn big(r: &Ratio<BigInt>) -> Qq {
    Qq::Big(r.clone())
}

impl Qq {
    /// Normalized rational, or `None` for a zero denominator.
    pub fn new(num: BigInt, den: BigInt) -> Option<Self> {
        if den.is_zero() {
            return None;
        }
        match (num.to_i64(), den.to_i64()) {
            (Some(n), Some(d)) => Some(small(n as i128, d as i128)),
            _ => Some(Qq::Big(Ratio::new(num, den))),
        }
    }

    /// Fast path for machine ints (the binding's hot seam).
    pub fn from_int_pair(num: i64, den: i64) -> Option<Self> {
        if den == 0 {
            return None;
        }
        Some(small(num as i128, den as i128))
    }

    /// Build from an `i64` without going through strings.
    pub fn from_int(n: i64) -> Self {
        Qq::Small(n, 1)
    }

    pub fn zero() -> Self {
        Qq::Small(0, 1)
    }

    pub fn one() -> Self {
        Qq::Small(1, 1)
    }

    pub fn is_zero(&self) -> bool {
        match self {
            Qq::Small(n, _) => *n == 0,
            Qq::Big(r) => r.is_zero(),
        }
    }

    pub fn is_one(&self) -> bool {
        match self {
            Qq::Small(n, d) => *n == 1 && *d == 1,
            Qq::Big(r) => r.is_one(),
        }
    }

    pub fn numer(&self) -> BigInt {
        match self {
            Qq::Small(n, _) => BigInt::from(*n),
            Qq::Big(r) => r.numer().clone(),
        }
    }

    pub fn denom(&self) -> BigInt {
        match self {
            Qq::Small(_, d) => BigInt::from(*d),
            Qq::Big(r) => r.denom().clone(),
        }
    }

    pub fn numer_string(&self) -> String {
        match self {
            Qq::Small(n, _) => n.to_string(),
            Qq::Big(r) => r.numer().to_string(),
        }
    }

    pub fn denom_string(&self) -> String {
        match self {
            Qq::Small(_, d) => d.to_string(),
            Qq::Big(r) => r.denom().to_string(),
        }
    }

    fn as_big(&self) -> Ratio<BigInt> {
        match self {
            Qq::Small(n, d) => Ratio::new(BigInt::from(*n), BigInt::from(*d)),
            Qq::Big(r) => r.clone(),
        }
    }

    fn add_small(an: i128, ad: i128, bn: i128, bd: i128, neg: bool) -> Qq {
        let bn = if neg { -bn } else { bn };
        // (an*bd + bn*ad) / (ad*bd), checked; overflow falls back to Big.
        let r = an
            .checked_mul(bd)
            .and_then(|x| bn.checked_mul(ad).and_then(|y| x.checked_add(y)))
            .and_then(|n| ad.checked_mul(bd).map(|d| (n, d)));
        match r {
            Some((n, d)) => small(n, d),
            None => {
                let (x, y) = (
                    Ratio::new(BigInt::from(an), BigInt::from(ad)),
                    Ratio::new(BigInt::from(bn), BigInt::from(bd)),
                );
                Qq::Big(if neg { x - y } else { x + y })
            }
        }
    }

    fn mul_small(an: i128, ad: i128, bn: i128, bd: i128) -> Qq {
        let r = an
            .checked_mul(bn)
            .and_then(|n| ad.checked_mul(bd).map(|d| (n, d)));
        match r {
            Some((n, d)) => small(n, d),
            None => Qq::Big(
                Ratio::new(BigInt::from(an), BigInt::from(ad))
                    * Ratio::new(BigInt::from(bn), BigInt::from(bd)),
            ),
        }
    }

    /// Exact quotient; `None` on division by zero.
    pub fn exact_div(&self, rhs: &Self) -> Option<Self> {
        match (self, rhs) {
            (Qq::Small(an, ad), Qq::Small(bn, bd)) => {
                if *bn == 0 {
                    return None;
                }
                // (an/bn...) = (an*bd)/(ad*bn), checked.
                let r = (*an as i128)
                    .checked_mul(*bd as i128)
                    .and_then(|n| (*ad as i128).checked_mul(*bn as i128).map(|d| (n, d)));
                match r {
                    Some((n, d)) => Some(small(n, d)),
                    None => Some(Qq::Big(self.as_big() / rhs.as_big())),
                }
            }
            _ => {
                if rhs.is_zero() {
                    return None;
                }
                Some(Qq::Big(self.as_big() / rhs.as_big()))
            }
        }
    }
}

impl std::ops::Add<&Qq> for &Qq {
    type Output = Qq;
    fn add(self, rhs: &Qq) -> Qq {
        match (self, rhs) {
            (Qq::Small(an, ad), Qq::Small(bn, bd)) => {
                Qq::add_small(*an as i128, *ad as i128, *bn as i128, *bd as i128, false)
            }
            _ => big(&(self.as_big() + rhs.as_big())),
        }
    }
}

impl std::ops::Sub<&Qq> for &Qq {
    type Output = Qq;
    fn sub(self, rhs: &Qq) -> Qq {
        match (self, rhs) {
            (Qq::Small(an, ad), Qq::Small(bn, bd)) => {
                Qq::add_small(*an as i128, *ad as i128, *bn as i128, *bd as i128, true)
            }
            _ => big(&(self.as_big() - rhs.as_big())),
        }
    }
}

impl std::ops::Mul<&Qq> for &Qq {
    type Output = Qq;
    fn mul(self, rhs: &Qq) -> Qq {
        match (self, rhs) {
            (Qq::Small(an, ad), Qq::Small(bn, bd)) => {
                Qq::mul_small(*an as i128, *ad as i128, *bn as i128, *bd as i128)
            }
            _ => big(&(self.as_big() * rhs.as_big())),
        }
    }
}

impl std::ops::Neg for &Qq {
    type Output = Qq;
    fn neg(self) -> Qq {
        match self {
            Qq::Small(n, d) => match n.checked_neg() {
                Some(m) => Qq::Small(m, *d),
                None => big(&-self.as_big()),
            },
            Qq::Big(_) => big(&-self.as_big()),
        }
    }
}

impl std::fmt::Display for Qq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Qq::Small(n, d) => {
                if *d == 1 {
                    write!(f, "{n}")
                } else {
                    write!(f, "{n}/{d}")
                }
            }
            Qq::Big(r) => {
                if r.is_integer() {
                    write!(f, "{}", r.numer())
                } else {
                    write!(f, "{}/{}", r.numer(), r.denom())
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn qq(n: i64, d: i64) -> Qq {
        Qq::from_int_pair(n, d).unwrap()
    }

    #[test]
    fn normalizes_sign_and_gcd() {
        assert_eq!(qq(2, 4), qq(1, 2));
        assert_eq!(qq(1, -2), qq(-1, 2));
        assert_eq!(qq(-3, -6), qq(1, 2));
        assert_eq!(qq(0, 7), Qq::zero());
    }

    #[test]
    fn rejects_zero_denominator() {
        assert!(Qq::new(1.into(), 0.into()).is_none());
        assert!(Qq::from_int_pair(1, 0).is_none());
    }

    #[test]
    fn arithmetic_matches_rationals() {
        assert_eq!(&qq(1, 2) + &qq(1, 3), qq(5, 6));
        assert_eq!(&qq(1, 2) - &qq(3, 4), qq(-1, 4));
        assert_eq!(&qq(2, 3) * &qq(9, 4), qq(3, 2));
        assert_eq!(qq(3, 4).exact_div(&qq(3, 2)), Some(qq(1, 2)));
        assert_eq!(qq(1, 2).exact_div(&Qq::zero()), None);
    }

    #[test]
    fn overflow_falls_back_to_big() {
        let big = &Qq::from_int(i64::MAX) + &Qq::from_int(1);
        assert!(matches!(big, Qq::Big(_)));
        assert_eq!(big.to_string(), "9223372036854775808");
        // and back down through division
        let back = big.exact_div(&Qq::from_int(2)).unwrap();
        assert_eq!(back.to_string(), "4611686018427387904");
    }

    #[test]
    fn big_values_round_trip() {
        let n = BigInt::parse_bytes(b"123456789012345678901234567890", 10).unwrap();
        let q = Qq::new(n.clone(), BigInt::from(7)).unwrap();
        assert!(matches!(q, Qq::Big(_)));
        // Ratio::new normalizes: n is divisible by 7.
        assert_eq!(q.numer(), n / 7);
        assert_eq!(q.denom(), BigInt::from(1));
    }
}
