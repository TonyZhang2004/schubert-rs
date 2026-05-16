use crate::{giambelli, pieri, Grassmannian, Partition, Result, SchubertError};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

/// A sparse integral linear combination of Schubert classes on one Grassmannian.
///
/// Expressions are stored in the Schubert basis as a deterministic
/// `BTreeMap<Partition, BigInt>`.  Arithmetic is checked against the embedded
/// [`Grassmannian`] context, so products such as `sigma1.pow(4)?` do not need a
/// separate ambient parameter.  Use [`Self::product`] for checked multiplication
/// and `*` when panicking on incompatible Grassmannians is acceptable.
///
/// ```
/// use num_bigint::BigInt;
/// use schubert_core::Grassmannian;
///
/// # fn main() -> schubert_core::Result<()> {
/// let g = Grassmannian::new(2, 4)?;
/// let h = g.class(vec![1])?;
/// let h2 = h.pow(2)?;
/// assert_eq!(h2.to_string(), "sigma_(1,1) + sigma_(2)");
/// assert_eq!(h.pow(4)?.integral(&g)?, BigInt::from(2));
/// # Ok(())
/// # }
/// ```
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SchubertExpr {
    g: Grassmannian,
    terms: BTreeMap<Partition, BigInt>,
}

impl SchubertExpr {
    /// The zero expression in `CH^*(g)`.
    pub fn zero(g: &Grassmannian) -> Self {
        Self {
            g: g.clone(),
            terms: BTreeMap::new(),
        }
    }

    /// The multiplicative unit `sigma_()` in `CH^*(g)`.
    pub fn one(g: &Grassmannian) -> Self {
        Self::from_class(g, Partition::zero(g))
    }

    /// Build a single Schubert basis class `sigma_partition` in `CH^*(g)`.
    ///
    /// The partition is revalidated against `g` when it appears in products,
    /// so prefer [`Grassmannian::class`] for user input.
    pub fn from_class(g: &Grassmannian, partition: Partition) -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(partition, BigInt::one());
        Self {
            g: g.clone(),
            terms,
        }
    }

    /// The ambient Grassmannian carried by this expression.
    pub fn grassmannian(&self) -> &Grassmannian {
        &self.g
    }

    /// Schubert-basis coefficients, ordered lexicographically by partition.
    pub fn terms(&self) -> &BTreeMap<Partition, BigInt> {
        &self.terms
    }

    /// Number of nonzero Schubert-basis terms.
    pub fn term_count(&self) -> usize {
        self.terms.len()
    }

    /// Whether all Schubert-basis coefficients are zero.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Whether this expression is the multiplicative unit `sigma_()`.
    pub fn is_one(&self) -> bool {
        self.terms.len() == 1
            && self
                .terms
                .iter()
                .next()
                .is_some_and(|(partition, coeff)| partition.is_empty() && coeff.is_one())
    }

    /// Coefficient of a Schubert basis class, or zero if it is absent.
    pub fn coefficient(&self, partition: &Partition) -> BigInt {
        self.terms
            .get(partition)
            .cloned()
            .unwrap_or_else(BigInt::zero)
    }

    /// Checked addition of two expressions on the same Grassmannian.
    pub fn checked_add(&self, rhs: &Self) -> Result<Self> {
        self.ensure_same_grassmannian(rhs)?;
        let mut out = self.clone();
        for (partition, coeff) in &rhs.terms {
            out.add_term(partition.clone(), coeff.clone());
        }
        Ok(out)
    }

    /// Checked subtraction of two expressions on the same Grassmannian.
    pub fn checked_sub(&self, rhs: &Self) -> Result<Self> {
        self.ensure_same_grassmannian(rhs)?;
        let mut out = self.clone();
        for (partition, coeff) in &rhs.terms {
            out.add_term(partition.clone(), -coeff);
        }
        Ok(out)
    }

    /// Multiply every Schubert coefficient by an integer scalar.
    pub fn scalar_mul(&self, scalar: &BigInt) -> Self {
        if scalar.is_zero() || self.is_zero() {
            return Self::zero(&self.g);
        }

        let terms = self
            .terms
            .iter()
            .map(|(partition, coeff)| (partition.clone(), coeff * scalar))
            .filter(|(_, coeff)| !coeff.is_zero())
            .collect();

        Self {
            g: self.g.clone(),
            terms,
        }
    }

    /// Multiply by the special Schubert class `sigma_r` using Pieri's rule.
    ///
    /// Out-of-range special classes (`r < 0` or `r > n-k`) are interpreted as
    /// zero, matching the Giambelli determinant convention.
    pub fn multiply_by_special(&self, r: isize) -> Result<Self> {
        if r < 0 || r as usize > self.g.width() {
            return Ok(Self::zero(&self.g));
        }

        if r == 0 {
            return Ok(self.clone());
        }

        let mut out = Self::zero(&self.g);
        for (lambda, coeff) in &self.terms {
            for nu in pieri::multiply_partition_by_special(&self.g, lambda, r as usize)? {
                out.add_term(nu, coeff.clone());
            }
        }
        Ok(out)
    }

    /// Multiply by one Schubert basis class `sigma_rhs` using Giambelli plus Pieri.
    pub fn mul_class(&self, rhs: &Partition) -> Result<Self> {
        rhs.validate_for(&self.g)?;
        giambelli::multiply_expr_by_class(self, rhs)
    }

    /// Schubert product of two expressions in the same Grassmannian.
    ///
    /// This is the clearest method name for public code.  The older [`Self::mul`]
    /// alias remains available for compatibility.
    pub fn product(&self, rhs: &Self) -> Result<Self> {
        self.ensure_same_grassmannian(rhs)?;

        let mut out = Self::zero(&self.g);
        for (partition, coeff) in &rhs.terms {
            let product = self.mul_class(partition)?.scalar_mul(coeff);
            out = out.checked_add(&product)?;
        }
        Ok(out)
    }

    /// Compatibility alias for [`Self::product`].
    pub fn mul(&self, rhs: &Self) -> Result<Self> {
        self.product(rhs)
    }

    /// Repeated Schubert product by exponentiation by squaring.
    pub fn pow(&self, exponent: u32) -> Result<Self> {
        let mut exp = exponent;
        let mut result = Self::one(&self.g);
        let mut base = self.clone();

        while exp > 0 {
            if exp % 2 == 1 {
                result = result.product(&base)?;
            }
            exp /= 2;
            if exp > 0 {
                base = base.product(&base)?;
            }
        }

        Ok(result)
    }

    /// Integrate over `g` by extracting the coefficient of the top class.
    ///
    /// The explicit `g` parameter is retained for API clarity and is checked
    /// against the expression's embedded context.
    pub fn integral(&self, g: &Grassmannian) -> Result<BigInt> {
        self.ensure_grassmannian(g)?;
        let top = g.top_partition()?;
        Ok(self.coefficient(&top))
    }

    pub(crate) fn ensure_grassmannian(&self, g: &Grassmannian) -> Result<()> {
        if &self.g == g {
            Ok(())
        } else {
            Err(SchubertError::IncompatibleGrassmannian {
                left_k: self.g.k(),
                left_n: self.g.n(),
                right_k: g.k(),
                right_n: g.n(),
            })
        }
    }

    pub(crate) fn add_term(&mut self, partition: Partition, coeff: BigInt) {
        if coeff.is_zero() {
            return;
        }

        let entry = self.terms.entry(partition).or_insert_with(BigInt::zero);
        *entry += coeff;
        if entry.is_zero() {
            self.terms.retain(|_, coeff| !coeff.is_zero());
        }
    }

    fn ensure_same_grassmannian(&self, rhs: &Self) -> Result<()> {
        if self.g == rhs.g {
            Ok(())
        } else {
            Err(SchubertError::IncompatibleGrassmannian {
                left_k: self.g.k(),
                left_n: self.g.n(),
                right_k: rhs.g.k(),
                right_n: rhs.g.n(),
            })
        }
    }
}

impl fmt::Display for SchubertExpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }

        for (index, (partition, coeff)) in self.terms.iter().enumerate() {
            let negative = coeff.is_negative();
            let abs_coeff = coeff.abs();

            if index == 0 {
                if negative {
                    write!(f, "-")?;
                }
            } else if negative {
                write!(f, " - ")?;
            } else {
                write!(f, " + ")?;
            }

            if partition.is_empty() {
                write!(f, "{abs_coeff}")?;
            } else if abs_coeff.is_one() {
                write!(f, "sigma_{partition}")?;
            } else {
                write!(f, "{abs_coeff}*sigma_{partition}")?;
            }
        }

        Ok(())
    }
}

impl Add for SchubertExpr {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.checked_add(&rhs)
            .expect("cannot add Schubert expressions on different Grassmannians")
    }
}

impl Sub for SchubertExpr {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.checked_sub(&rhs)
            .expect("cannot subtract Schubert expressions on different Grassmannians")
    }
}

impl Neg for SchubertExpr {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.scalar_mul(&BigInt::from(-1))
    }
}

impl Mul for SchubertExpr {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.product(&rhs)
            .expect("cannot multiply Schubert expressions on different Grassmannians")
    }
}

impl std::ops::Mul<BigInt> for SchubertExpr {
    type Output = Self;

    fn mul(self, rhs: BigInt) -> Self::Output {
        self.scalar_mul(&rhs)
    }
}

impl std::ops::Mul<i64> for SchubertExpr {
    type Output = Self;

    fn mul(self, rhs: i64) -> Self::Output {
        self.scalar_mul(&BigInt::from(rhs))
    }
}
