//! Python bindings for `sympy-poly-core` (PyO3).
//!
//! Seam convention: coefficients cross as `(num, den)` int pairs in
//! sympy's high-to-low dense order; the binding converts to the core's
//! low-to-high order and back. Rationals return as `(num_str, den_str)`
//! string pairs (exact; the `sympy_poly` wrapper turns them into `QQ`).

use num_bigint::BigInt;
use pyo3::exceptions::{PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use std::str::FromStr;
use sympy_poly_core::{Dmp, Qq};

/// Python int -> `Qq` fast path: machine ints stay machine ints; only
/// genuinely big values pay for decimal spelling + parse.
fn to_bigint(obj: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    if let Ok(n) = obj.extract::<i64>() {
        return Ok(BigInt::from(n));
    }
    let s: String = obj.str()?.extract()?;
    BigInt::from_str(&s).map_err(|e| PyValueError::new_err(format!("not an int: {e}")))
}

fn to_qq(pair: &Bound<'_, PyAny>) -> PyResult<Qq> {
    let (n, d): (Bound<'_, PyAny>, Bound<'_, PyAny>) = pair.extract().map_err(|_| {
        PyValueError::new_err("coefficient must be a (num, den) pair of ints")
    })?;
    if let (Ok(n), Ok(d)) = (n.extract::<i64>(), d.extract::<i64>()) {
        return Qq::from_int_pair(n, d)
            .ok_or_else(|| PyZeroDivisionError::new_err("zero denominator in QQ"));
    }
    let (n, d) = (to_bigint(&n)?, to_bigint(&d)?);
    Qq::new(n, d).ok_or_else(|| PyZeroDivisionError::new_err("zero denominator in QQ"))
}

fn qq_pair(q: &Qq) -> (String, String) {
    (q.numer_string(), q.denom_string())
}

/// High-to-low pairs in -> low-to-high core vec.
fn to_dup(items: Vec<Bound<'_, PyAny>>) -> PyResult<Vec<Qq>> {
    let mut v: Vec<Qq> = items.iter().map(to_qq).collect::<PyResult<_>>()?;
    v.reverse();
    while v.last().is_some_and(Qq::is_zero) {
        v.pop();
    }
    Ok(v)
}

/// Low-to-high core vec -> high-to-low string pairs out.
fn from_dup(mut v: Vec<Qq>) -> Vec<(String, String)> {
    while v.last().is_some_and(Qq::is_zero) {
        v.pop();
    }
    v.iter().rev().map(qq_pair).collect()
}

/// Nested high-to-low lists in -> core `Dmp` at `depth` levels.
fn to_dmp(obj: &Bound<'_, PyAny>, depth: u32) -> PyResult<Dmp> {
    if depth == 0 {
        return Ok(Dmp::Coeff(to_qq(obj)?));
    }
    let items: Vec<Bound<'_, PyAny>> = obj.extract().map_err(|_| {
        PyValueError::new_err("dmp level must be a list (use [] for the zero polynomial)")
    })?;
    let mut cs: Vec<Dmp> = items
        .iter()
        .map(|c| to_dmp(c, depth - 1))
        .collect::<PyResult<_>>()?;
    cs.reverse();
    Ok(Dmp::Poly(cs).strip())
}

fn from_dmp(d: &Dmp, out: &mut Vec<PyObject>, py: Python<'_>) {
    match d {
        Dmp::Coeff(q) => {
            let (n, s) = qq_pair(q);
            out.push((n, s).into_pyobject(py).unwrap().into());
        }
        Dmp::Poly(cs) => {
            let mut level = Vec::with_capacity(cs.len());
            for c in cs.iter().rev() {
                let mut child = Vec::new();
                from_dmp(c, &mut child, py);
                debug_assert_eq!(child.len(), 1);
                level.push(child.pop().unwrap());
            }
            out.push(level.into_pyobject(py).unwrap().into());
        }
    }
}

fn export_dmp(py: Python<'_>, d: &Dmp) -> PyResult<PyObject> {
    let mut out = Vec::new();
    from_dmp(d, &mut out, py);
    debug_assert_eq!(out.len(), 1);
    Ok(out.pop().unwrap())
}

#[pyfunction]
fn dup_mul(f: Vec<Bound<'_, PyAny>>, g: Vec<Bound<'_, PyAny>>) -> PyResult<Vec<(String, String)>> {
    Ok(from_dup(sympy_poly_core::dup_mul(&to_dup(f)?, &to_dup(g)?)))
}

#[pyfunction]
fn dup_add(f: Vec<Bound<'_, PyAny>>, g: Vec<Bound<'_, PyAny>>) -> PyResult<Vec<(String, String)>> {
    Ok(from_dup(sympy_poly_core::dup_add(&to_dup(f)?, &to_dup(g)?)))
}

#[pyfunction]
fn dup_sub(f: Vec<Bound<'_, PyAny>>, g: Vec<Bound<'_, PyAny>>) -> PyResult<Vec<(String, String)>> {
    Ok(from_dup(sympy_poly_core::dup_sub(&to_dup(f)?, &to_dup(g)?)))
}

#[pyfunction]
fn dup_gcd(f: Vec<Bound<'_, PyAny>>, g: Vec<Bound<'_, PyAny>>) -> PyResult<Vec<(String, String)>> {
    Ok(from_dup(sympy_poly_core::dup_gcd(&to_dup(f)?, &to_dup(g)?)))
}

#[pyfunction]
fn dmp_mul(a: Bound<'_, PyAny>, b: Bound<'_, PyAny>, levels: u32) -> PyResult<PyObject> {
    let depth = levels;
    if levels == 0 {
        return Err(PyValueError::new_err("dmp_mul needs levels >= 1"));
    }
    let h = sympy_poly_core::dmp_mul(&to_dmp(&a, depth)?, &to_dmp(&b, depth)?, depth);
    export_dmp(a.py(), &h)
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dup_mul, m)?)?;
    m.add_function(wrap_pyfunction!(dup_add, m)?)?;
    m.add_function(wrap_pyfunction!(dup_sub, m)?)?;
    m.add_function(wrap_pyfunction!(dup_gcd, m)?)?;
    m.add_function(wrap_pyfunction!(dmp_mul, m)?)?;
    Ok(())
}
