use crate::{Partition, Result, SchubertError, SchubertExpr};
use serde::{Deserialize, Serialize};

/// Shape data for the ordinary Grassmannian `Gr(k,n)`.
///
/// The Schubert basis of `CH^*(Gr(k,n), Z)` is indexed by partitions fitting
/// inside the `k x (n-k)` rectangle.  The stored `width` is `n-k`, so the top
/// Schubert class is indexed by `(n-k, ..., n-k)`.
///
/// ```
/// use schubert_core::Grassmannian;
///
/// # fn main() -> schubert_core::Result<()> {
/// let g = Grassmannian::new(2, 5)?;
/// assert_eq!(g.k(), 2);
/// assert_eq!(g.width(), 3);
/// assert_eq!(g.dimension(), 6);
/// assert_eq!(g.top_partition()?.parts(), &[3, 3]);
/// # Ok(())
/// # }
/// ```
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Grassmannian {
    k: usize,
    n: usize,
    w: usize,
}

impl Grassmannian {
    /// Construct `Gr(k,n)`.
    ///
    /// Returns an error unless `0 < k <= n`.  The rectangle width `n-k` must
    /// also fit in `u16`, matching the exact storage used by [`Partition`].
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

    /// The dimension of the chosen subspaces in `Gr(k,n)`.
    pub fn k(&self) -> usize {
        self.k
    }

    /// The ambient vector-space dimension in `Gr(k,n)`.
    pub fn n(&self) -> usize {
        self.n
    }

    /// The width `n-k` of the indexing rectangle.
    pub fn width(&self) -> usize {
        self.w
    }

    /// The complex dimension `k(n-k)` of the Grassmannian.
    ///
    /// This is also the codimension of the top Schubert class.
    pub fn dimension(&self) -> usize {
        self.k * self.w
    }

    /// The rectangular partition `(n-k, ..., n-k)` indexing the point class.
    pub fn top_partition(&self) -> Result<Partition> {
        Partition::new(vec![self.w as u16; self.k], self)
    }

    /// Build the Schubert basis element `sigma_lambda` for `lambda = parts`.
    ///
    /// Input is normalized by padding trailing zeroes to length `k`.  For
    /// example, in `Gr(2,4)`, `vec![1]` and `vec![1, 0]` both denote
    /// `sigma_(1)`.
    ///
    /// ```
    /// use schubert_core::Grassmannian;
    ///
    /// # fn main() -> schubert_core::Result<()> {
    /// let g = Grassmannian::new(2, 4)?;
    /// let sigma_1 = g.class(vec![1])?;
    /// assert_eq!(sigma_1.to_string(), "sigma_(1)");
    /// # Ok(())
    /// # }
    /// ```
    pub fn class(&self, parts: Vec<u16>) -> Result<SchubertExpr> {
        let partition = Partition::new(parts, self)?;
        Ok(SchubertExpr::from_class(self, partition))
    }
}
