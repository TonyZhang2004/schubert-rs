use thiserror::Error;

/// Errors returned by checked Schubert-calculus operations.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SchubertError {
    /// `Gr(k,n)` requires `0 < k <= n`.
    #[error("invalid Grassmannian Gr({k},{n}); expected 0 < k <= n")]
    InvalidGrassmannian { k: usize, n: usize },

    /// The rectangle width `n-k` cannot be represented by `u16` partition parts.
    #[error("Grassmannian width {width} is too large for u16 partition parts")]
    WidthTooLarge { width: usize },

    /// A partition has more than `k` rows.
    #[error("partition has length {len}, but Gr({k},{n}) allows at most {k} parts")]
    PartitionTooLong { len: usize, k: usize, n: usize },

    /// A partition part lies outside the `k x (n-k)` rectangle.
    #[error("partition part {part} exceeds rectangle width {width} for Gr({k},{n})")]
    PartitionPartTooLarge {
        part: u16,
        width: usize,
        k: usize,
        n: usize,
    },

    /// A partition is not weakly decreasing.
    #[error("partition is not weakly decreasing at index {index}: {left} < {right}")]
    PartitionNotWeaklyDecreasing { index: usize, left: u16, right: u16 },

    /// An operation mixed expressions from different Grassmannians.
    #[error(
        "incompatible Grassmannians: left Gr({left_k},{left_n}), right Gr({right_k},{right_n})"
    )]
    IncompatibleGrassmannian {
        left_k: usize,
        left_n: usize,
        right_k: usize,
        right_n: usize,
    },

    /// A normalized partition has the wrong length for the target Grassmannian.
    #[error("partition length {len} is incompatible with Gr({k},{n})")]
    PartitionWrongLength { len: usize, k: usize, n: usize },
}
