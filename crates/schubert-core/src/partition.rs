use crate::{Grassmannian, Result, SchubertError};
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Partition {
    parts: Vec<u16>,
}

impl Partition {
    pub fn new(parts: Vec<u16>, g: &Grassmannian) -> Result<Self> {
        Self::from_parts(parts, g.k(), g.width(), g.k(), g.n())
    }

    pub(crate) fn zero(g: &Grassmannian) -> Self {
        Self {
            parts: vec![0; g.k()],
        }
    }

    pub(crate) fn from_parts(
        mut parts: Vec<u16>,
        k: usize,
        w: usize,
        g_k: usize,
        g_n: usize,
    ) -> Result<Self> {
        if parts.len() > k {
            return Err(SchubertError::PartitionTooLong {
                len: parts.len(),
                k: g_k,
                n: g_n,
            });
        }

        for i in 1..parts.len() {
            if parts[i - 1] < parts[i] {
                return Err(SchubertError::PartitionNotWeaklyDecreasing {
                    index: i - 1,
                    left: parts[i - 1],
                    right: parts[i],
                });
            }
        }

        for &part in &parts {
            if part as usize > w {
                return Err(SchubertError::PartitionPartTooLarge {
                    part,
                    width: w,
                    k: g_k,
                    n: g_n,
                });
            }
        }

        parts.resize(k, 0);
        Ok(Self { parts })
    }

    pub(crate) fn from_normalized_unchecked(parts: Vec<u16>) -> Self {
        Self { parts }
    }

    pub fn parts(&self) -> &[u16] {
        &self.parts
    }

    pub fn size(&self) -> usize {
        self.parts.iter().map(|&part| part as usize).sum()
    }

    pub fn is_top(&self, g: &Grassmannian) -> bool {
        self.parts.len() == g.k() && self.parts.iter().all(|&part| part as usize == g.width())
    }

    pub(crate) fn validate_for(&self, g: &Grassmannian) -> Result<()> {
        if self.parts.len() != g.k() {
            return Err(SchubertError::PartitionWrongLength {
                len: self.parts.len(),
                k: g.k(),
                n: g.n(),
            });
        }

        for i in 1..self.parts.len() {
            if self.parts[i - 1] < self.parts[i] {
                return Err(SchubertError::PartitionNotWeaklyDecreasing {
                    index: i - 1,
                    left: self.parts[i - 1],
                    right: self.parts[i],
                });
            }
        }

        for &part in &self.parts {
            if part as usize > g.width() {
                return Err(SchubertError::PartitionPartTooLarge {
                    part,
                    width: g.width(),
                    k: g.k(),
                    n: g.n(),
                });
            }
        }

        Ok(())
    }
}

impl fmt::Display for Partition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let visible_len = self
            .parts
            .iter()
            .rposition(|&part| part != 0)
            .map_or(0, |index| index + 1);

        write!(f, "(")?;
        for (index, part) in self.parts[..visible_len].iter().enumerate() {
            if index > 0 {
                write!(f, ",")?;
            }
            write!(f, "{part}")?;
        }
        write!(f, ")")
    }
}
