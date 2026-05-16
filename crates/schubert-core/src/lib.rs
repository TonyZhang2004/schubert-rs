//! Exact Schubert calculus for ordinary integral cohomology of type-A Grassmannians.
//!
//! This crate computes in the Schubert basis of `CH^*(Gr(k,n), Z)`.  A
//! [`Grassmannian`] fixes the ambient `k x (n-k)` rectangle, [`Partition`] values
//! index Schubert classes `sigma_lambda`, and [`SchubertExpr`] stores sparse
//! integral Schubert expansions.  Multiplication uses direct
//! Littlewood-Richardson tableau counting for general products and Pieri's rule
//! for special classes.
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
//! assert_eq!(ans.to_string(), "2*sigma_(2,2)");
//! # Ok(())
//! # }
//! ```

mod errors;
mod expr;
#[cfg(test)]
mod giambelli;
mod grassmannian;
mod littlewood_richardson;
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

    fn partitions(g: &Grassmannian) -> Vec<Partition> {
        fn rec(
            g: &Grassmannian,
            row: usize,
            max_part: u16,
            current: &mut Vec<u16>,
            out: &mut Vec<Partition>,
        ) {
            if row == g.k() {
                out.push(Partition::new(current.clone(), g).unwrap());
                return;
            }

            for part in (0..=max_part).rev() {
                current.push(part);
                rec(g, row + 1, part, current, out);
                current.pop();
            }
        }

        let mut out = Vec::new();
        rec(g, 0, g.width() as u16, &mut Vec::new(), &mut out);
        out
    }

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
    fn expression_display_and_product_alias_are_ergonomic() {
        let g = Grassmannian::new(2, 4).unwrap();
        let sigma_1 = g.class(vec![1]).unwrap();
        let sigma_2 = g.class(vec![2]).unwrap();

        let square = sigma_1.product(&sigma_1).unwrap();
        assert_eq!(square.to_string(), "sigma_(1,1) + sigma_(2)");
        assert_eq!(square.term_count(), 2);
        assert_eq!(sigma_1.clone() * sigma_1, square);

        let mixed = square.checked_sub(&sigma_2).unwrap();
        assert_eq!(mixed.to_string(), "sigma_(1,1)");
        assert!(!mixed.is_one());
        assert_eq!(SchubertExpr::one(&g).to_string(), "1");
        assert!(SchubertExpr::one(&g).is_one());
        assert_eq!(SchubertExpr::zero(&g).to_string(), "0");
    }

    #[test]
    fn projective_space_relation_for_gr_1_n() {
        let g = Grassmannian::new(1, 5).unwrap();
        let h = g.class(vec![1]).unwrap();

        assert_eq!(h.pow(4).unwrap().integral(&g).unwrap(), BigInt::one());
        assert_eq!(h.pow(5).unwrap().integral(&g).unwrap(), BigInt::zero());
        assert!(h.pow(5).unwrap().is_zero());
    }

    #[test]
    fn littlewood_richardson_matches_giambelli_in_gr_2_4() {
        let g = Grassmannian::new(2, 4).unwrap();
        let all = partitions(&g);

        for lambda in &all {
            let expr = SchubertExpr::from_class(&g, lambda.clone());
            for mu in &all {
                let lr = expr.mul_class(mu).unwrap();
                let giambelli = giambelli::multiply_expr_by_class(&expr, mu).unwrap();
                assert_eq!(lr, giambelli, "lambda={lambda}, mu={mu}");
            }
        }
    }

    #[test]
    fn littlewood_richardson_matches_giambelli_in_gr_3_5() {
        let g = Grassmannian::new(3, 5).unwrap();
        let all = partitions(&g);

        for lambda in &all {
            let expr = SchubertExpr::from_class(&g, lambda.clone());
            for mu in &all {
                let lr = expr.mul_class(mu).unwrap();
                let giambelli = giambelli::multiply_expr_by_class(&expr, mu).unwrap();
                assert_eq!(lr, giambelli, "lambda={lambda}, mu={mu}");
            }
        }
    }

    #[test]
    fn littlewood_richardson_special_classes_match_pieri() {
        let g = Grassmannian::new(3, 5).unwrap();

        for lambda in partitions(&g) {
            let expr = SchubertExpr::from_class(&g, lambda);
            for r in 0..=g.width() {
                let special = Partition::new(vec![r as u16], &g).unwrap();
                let lr = expr.mul_class(&special).unwrap();
                let pieri = expr.multiply_by_special(r as isize).unwrap();
                assert_eq!(lr, pieri);
            }
        }
    }

    #[test]
    fn littlewood_richardson_unit_and_overfull_products() {
        let g = Grassmannian::new(2, 4).unwrap();
        let sigma_2 = g.class(vec![2]).unwrap();
        let one = Partition::new(vec![], &g).unwrap();
        let top = g.top_partition().unwrap();

        assert_eq!(sigma_2.mul_class(&one).unwrap(), sigma_2);
        assert!(sigma_2.mul_class(&top).unwrap().is_zero());
    }

    #[test]
    fn littlewood_richardson_counts_coefficients_above_one() {
        let g = Grassmannian::new(3, 6).unwrap();
        let lambda = Partition::new(vec![2, 1], &g).unwrap();
        let mu = Partition::new(vec![2, 1], &g).unwrap();
        let nu = Partition::new(vec![3, 2, 1], &g).unwrap();

        let product = SchubertExpr::from_class(&g, lambda).mul_class(&mu).unwrap();

        assert_eq!(product.coefficient(&nu), BigInt::from(2));
    }

    #[test]
    fn littlewood_richardson_full_square_of_21_in_large_rectangle() {
        let g = Grassmannian::new(4, 8).unwrap();
        let sigma_21 = g.class(vec![2, 1]).unwrap();
        let product = sigma_21.product(&sigma_21).unwrap();

        let expected = [
            (vec![4, 2], 1),
            (vec![4, 1, 1], 1),
            (vec![3, 3], 1),
            (vec![3, 2, 1], 2),
            (vec![3, 1, 1, 1], 1),
            (vec![2, 2, 2], 1),
            (vec![2, 2, 1, 1], 1),
        ];

        assert_eq!(product.term_count(), expected.len());
        for (parts, coefficient) in expected {
            let partition = Partition::new(parts, &g).unwrap();
            assert_eq!(product.coefficient(&partition), BigInt::from(coefficient));
        }
    }

    #[test]
    fn littlewood_richardson_square_of_21_truncates_to_gr_3_6_rectangle() {
        let g = Grassmannian::new(3, 6).unwrap();
        let sigma_21 = g.class(vec![2, 1]).unwrap();
        let product = sigma_21.product(&sigma_21).unwrap();

        let expected = [(vec![3, 3], 1), (vec![3, 2, 1], 2), (vec![2, 2, 2], 1)];

        assert_eq!(product.term_count(), expected.len());
        for (parts, coefficient) in expected {
            let partition = Partition::new(parts, &g).unwrap();
            assert_eq!(product.coefficient(&partition), BigInt::from(coefficient));
        }
    }

    #[test]
    fn littlewood_richardson_mixed_product_has_four_hand_checked_outputs() {
        let g = Grassmannian::new(3, 7).unwrap();
        let sigma_22 = g.class(vec![2, 2]).unwrap();
        let sigma_21 = g.class(vec![2, 1]).unwrap();
        let product = sigma_22.product(&sigma_21).unwrap();

        let expected = [vec![4, 3], vec![4, 2, 1], vec![3, 3, 1], vec![3, 2, 2]];

        assert_eq!(product.term_count(), expected.len());
        for parts in expected {
            let partition = Partition::new(parts, &g).unwrap();
            assert_eq!(product.coefficient(&partition), BigInt::one());
        }
    }
}
