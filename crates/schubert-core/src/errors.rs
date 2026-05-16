use thiserror::Error;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SchubertError {
    #[error("invalid Grassmannian Gr({k},{n}); expected 0 < k <= n")]
    InvalidGrassmannian { k: usize, n: usize },

    #[error("Grassmannian width {width} is too large for u16 partition parts")]
    WidthTooLarge { width: usize },

    #[error("partition has length {len}, but Gr({k},{n}) allows at most {k} parts")]
    PartitionTooLong { len: usize, k: usize, n: usize },

    #[error("partition part {part} exceeds rectangle width {width} for Gr({k},{n})")]
    PartitionPartTooLarge {
        part: u16,
        width: usize,
        k: usize,
        n: usize,
    },

    #[error("partition is not weakly decreasing at index {index}: {left} < {right}")]
    PartitionNotWeaklyDecreasing { index: usize, left: u16, right: u16 },

    #[error(
        "incompatible Grassmannians: left Gr({left_k},{left_n}), right Gr({right_k},{right_n})"
    )]
    IncompatibleGrassmannian {
        left_k: usize,
        left_n: usize,
        right_k: usize,
        right_n: usize,
    },

    #[error("partition length {len} is incompatible with Gr({k},{n})")]
    PartitionWrongLength { len: usize, k: usize, n: usize },
}
