use crate::{Partition, Result, SchubertError, SchubertExpr};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Grassmannian {
    k: usize,
    n: usize,
    w: usize,
}

impl Grassmannian {
    pub fn new(k: usize, n: usize) -> Result<Self> {
        if k == 0 || k > n {
            return Err(SchubertError::InvalidGrassmannian { k, n });
        }

        let w = n - k;
        if w > u16::MAX as usize {
            return Err(SchubertError::WidthTooLarge { width: w });
        }

        Ok(Self { k, n, w })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn width(&self) -> usize {
        self.w
    }

    pub fn dimension(&self) -> usize {
        self.k * self.w
    }

    pub fn top_partition(&self) -> Result<Partition> {
        Partition::new(vec![self.w as u16; self.k], self)
    }

    pub fn class(&self, parts: Vec<u16>) -> Result<SchubertExpr> {
        let partition = Partition::new(parts, self)?;
        Ok(SchubertExpr::from_class(self, partition))
    }
}
