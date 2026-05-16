use crate::{giambelli, pieri, Grassmannian, Partition, Result, SchubertError};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::ops::{Add, Neg, Sub};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SchubertExpr {
    g: Grassmannian,
    terms: BTreeMap<Partition, BigInt>,
}

impl SchubertExpr {
    pub fn zero(g: &Grassmannian) -> Self {
        Self {
            g: g.clone(),
            terms: BTreeMap::new(),
        }
    }

    pub fn one(g: &Grassmannian) -> Self {
        Self::from_class(g, Partition::zero(g))
    }

    pub fn from_class(g: &Grassmannian, partition: Partition) -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(partition, BigInt::one());
        Self {
            g: g.clone(),
            terms,
        }
    }

    pub fn grassmannian(&self) -> &Grassmannian {
        &self.g
    }

    pub fn terms(&self) -> &BTreeMap<Partition, BigInt> {
        &self.terms
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn coefficient(&self, partition: &Partition) -> BigInt {
        self.terms
            .get(partition)
            .cloned()
            .unwrap_or_else(BigInt::zero)
    }

    pub fn checked_add(&self, rhs: &Self) -> Result<Self> {
        self.ensure_same_grassmannian(rhs)?;
        let mut out = self.clone();
        for (partition, coeff) in &rhs.terms {
            out.add_term(partition.clone(), coeff.clone());
        }
        Ok(out)
    }

    pub fn checked_sub(&self, rhs: &Self) -> Result<Self> {
        self.ensure_same_grassmannian(rhs)?;
        let mut out = self.clone();
        for (partition, coeff) in &rhs.terms {
            out.add_term(partition.clone(), -coeff);
        }
        Ok(out)
    }

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

    pub fn mul_class(&self, rhs: &Partition) -> Result<Self> {
        rhs.validate_for(&self.g)?;
        giambelli::multiply_expr_by_class(self, rhs)
    }

    pub fn mul(&self, rhs: &Self) -> Result<Self> {
        self.ensure_same_grassmannian(rhs)?;

        let mut out = Self::zero(&self.g);
        for (partition, coeff) in &rhs.terms {
            let product = self.mul_class(partition)?.scalar_mul(coeff);
            out = out.checked_add(&product)?;
        }
        Ok(out)
    }

    pub fn pow(&self, exponent: u32) -> Result<Self> {
        let mut exp = exponent;
        let mut result = Self::one(&self.g);
        let mut base = self.clone();

        while exp > 0 {
            if exp % 2 == 1 {
                result = result.mul(&base)?;
            }
            exp /= 2;
            if exp > 0 {
                base = base.mul(&base)?;
            }
        }

        Ok(result)
    }

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
