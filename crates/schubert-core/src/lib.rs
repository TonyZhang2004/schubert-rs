//! Exact Schubert calculus for ordinary integral cohomology of type-A Grassmannians.
//!
//! ```
//! use num_bigint::BigInt;
//! use schubert_core::Grassmannian;
//!
//! # fn main() -> schubert_core::Result<()> {
//! let g = Grassmannian::new(2, 4)?;
//! let sigma1 = g.class(vec![1])?;
//! let ans = sigma1.pow(4)?;
//! assert_eq!(ans.integral(&g)?, BigInt::from(2));
//! # Ok(())
//! # }
//! ```

mod errors;
mod expr;
mod giambelli;
mod grassmannian;
mod partition;
mod pieri;

pub use errors::SchubertError;
pub use expr::SchubertExpr;
pub use grassmannian::Grassmannian;
pub use partition::Partition;

pub type Result<T> = std::result::Result<T, SchubertError>;

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_traits::{One, Zero};

    #[test]
    fn grassmannian_validation_and_helpers() {
        assert!(Grassmannian::new(0, 4).is_err());
        assert!(Grassmannian::new(5, 4).is_err());

        let g = Grassmannian::new(2, 4).unwrap();
        assert_eq!(g.k(), 2);
        assert_eq!(g.n(), 4);
        assert_eq!(g.width(), 2);
        assert_eq!(g.dimension(), 4);
        assert_eq!(g.top_partition().unwrap().parts(), &[2, 2]);
    }

    #[test]
    fn partition_normalization_and_validation() {
        let g = Grassmannian::new(3, 6).unwrap();
        let p = Partition::new(vec![2, 1], &g).unwrap();
        assert_eq!(p.parts(), &[2, 1, 0]);
        assert_eq!(p.size(), 3);

        assert!(Partition::new(vec![1, 2], &g).is_err());
        assert!(Partition::new(vec![4], &g).is_err());
        assert!(Partition::new(vec![1, 1, 1, 1], &g).is_err());
    }

    #[test]
    fn pieri_identity_and_out_of_range_special_class() {
        let g = Grassmannian::new(2, 4).unwrap();
        let sigma_11 = g.class(vec![1, 1]).unwrap();

        assert_eq!(sigma_11.multiply_by_special(0).unwrap(), sigma_11);
        assert!(sigma_11.multiply_by_special(3).unwrap().is_zero());
        assert!(sigma_11.multiply_by_special(-1).unwrap().is_zero());
    }

    #[test]
    fn pieri_basic_products_in_gr_2_4() {
        let g = Grassmannian::new(2, 4).unwrap();
        let sigma_1 = g.class(vec![1]).unwrap();

        let square = sigma_1.multiply_by_special(1).unwrap();
        assert_eq!(
            square
                .coefficient(&Partition::new(vec![2], &g).unwrap())
                .clone(),
            BigInt::one()
        );
        assert_eq!(
            square
                .coefficient(&Partition::new(vec![1, 1], &g).unwrap())
                .clone(),
            BigInt::one()
        );

        let sigma_2 = g.class(vec![2]).unwrap();
        let top_product = sigma_2.multiply_by_special(2).unwrap();
        assert_eq!(
            top_product.coefficient(&g.top_partition().unwrap()).clone(),
            BigInt::one()
        );
    }

    #[test]
    fn giambelli_multiplies_non_special_class() {
        let g = Grassmannian::new(2, 4).unwrap();
        let sigma_1 = g.class(vec![1]).unwrap();
        let sigma_11 = Partition::new(vec![1, 1], &g).unwrap();

        let product = sigma_1.mul_class(&sigma_11).unwrap();
        assert_eq!(
            product
                .coefficient(&Partition::new(vec![2, 1], &g).unwrap())
                .clone(),
            BigInt::one()
        );
        assert_eq!(product.terms().len(), 1);
    }

    #[test]
    fn required_gr_2_4_example() {
        let g = Grassmannian::new(2, 4).unwrap();
        let sigma1 = g.class(vec![1]).unwrap();
        let ans = sigma1.pow(4).unwrap();

        assert_eq!(ans.integral(&g).unwrap(), BigInt::from(2));
    }

    #[test]
    fn small_commutativity_and_associativity_checks() {
        let g = Grassmannian::new(2, 4).unwrap();
        let a = g.class(vec![1]).unwrap();
        let b = g.class(vec![2]).unwrap();
        let c = g.class(vec![1, 1]).unwrap();

        assert_eq!(a.mul(&b).unwrap(), b.mul(&a).unwrap());

        let left = a.mul(&b).unwrap().mul(&c).unwrap();
        let right = a.mul(&b.mul(&c).unwrap()).unwrap();
        assert_eq!(left, right);
    }

    #[test]
    fn projective_space_relation_for_gr_1_n() {
        let g = Grassmannian::new(1, 5).unwrap();
        let h = g.class(vec![1]).unwrap();

        assert_eq!(h.pow(4).unwrap().integral(&g).unwrap(), BigInt::one());
        assert_eq!(h.pow(5).unwrap().integral(&g).unwrap(), BigInt::zero());
        assert!(h.pow(5).unwrap().is_zero());
    }
}
