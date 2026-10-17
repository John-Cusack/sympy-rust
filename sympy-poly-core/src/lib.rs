//! `sympy-poly-core`: dense polynomial arithmetic over `QQ`.
//!
//! Reusable Rust core with no Python dependency. It ports the proven
//! seam of sympy's polys domains (headroom: polys-domains 0.85 rewrite):
//! the `QQ` domain operations and the dense `dup_*` / recursive `dmp_*`
//! arithmetic (schoolbook multiply, exact division, monic GCD).
pub mod dense;
pub mod qq;
pub mod zz;

pub use dense::{Dmp, dmp_add, dmp_mul, dup_add, dup_div, dup_gcd, dup_monic, dup_mul, dup_neg, dup_sub};
pub use qq::Qq;
