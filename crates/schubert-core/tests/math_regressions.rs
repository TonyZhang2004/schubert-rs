use num_bigint::BigInt;
use num_traits::{One, Zero};
use schubert_core::{Grassmannian, Partition, SchubertExpr};

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

fn dual_partition(lambda: &Partition, g: &Grassmannian) -> Partition {
    let w = g.width() as u16;
    let parts = lambda
        .parts()
        .iter()
        .rev()
        .map(|part| w - part)
        .collect::<Vec<_>>();
    Partition::new(parts, g).unwrap()
}

fn class(g: &Grassmannian, parts: &[u16]) -> SchubertExpr {
    g.class(parts.to_vec()).unwrap()
}

#[test]
fn gr_1_n_has_truncated_polynomial_behavior() {
    // Gr(1,n) is projective space, so its Schubert ring is Z[h]/(h^n).
    let g = Grassmannian::new(1, 6).unwrap();
    let h = class(&g, &[1]);

    for r in 0..=5 {
        let expected = Partition::new(vec![r], &g).unwrap();
        let power = h.pow(r as u32).unwrap();
        assert_eq!(power.terms().len(), 1);
        assert_eq!(power.coefficient(&expected), BigInt::one());
    }

    assert!(h.pow(6).unwrap().is_zero());
}

#[test]
fn gr_2_4_has_expected_low_degree_products() {
    // These are the Pieri paths inside the 2 x 2 rectangle for Gr(2,4).
    let g = Grassmannian::new(2, 4).unwrap();
    let h = class(&g, &[1]);
    let sigma_2 = Partition::new(vec![2], &g).unwrap();
    let sigma_11 = Partition::new(vec![1, 1], &g).unwrap();
    let sigma_21 = Partition::new(vec![2, 1], &g).unwrap();
    let top = g.top_partition().unwrap();

    let h2 = h.pow(2).unwrap();
    assert_eq!(h2.coefficient(&sigma_2), BigInt::one());
    assert_eq!(h2.coefficient(&sigma_11), BigInt::one());
    assert_eq!(h2.terms().len(), 2);

    let h3 = h.pow(3).unwrap();
    assert_eq!(h3.coefficient(&sigma_21), BigInt::from(2));
    assert_eq!(h3.terms().len(), 1);

    let h4 = h.pow(4).unwrap();
    assert_eq!(h4.coefficient(&top), BigInt::from(2));
    assert_eq!(h4.terms().len(), 1);

    assert_eq!(
        class(&g, &[2]).pow(2).unwrap().coefficient(&top),
        BigInt::one()
    );
    assert_eq!(
        class(&g, &[1, 1]).pow(2).unwrap().coefficient(&top),
        BigInt::one()
    );
    assert!(class(&g, &[2]).mul(&class(&g, &[1, 1])).unwrap().is_zero());
}

#[test]
fn poincare_duality_pairing_is_kronecker_delta_in_gr_2_4() {
    // Integration gives the Schubert Poincare pairing against the rotated complement.
    let g = Grassmannian::new(2, 4).unwrap();
    let all = partitions(&g);

    for lambda in &all {
        let dual = dual_partition(lambda, &g);
        for mu in &all {
            let expected = if mu == &dual {
                BigInt::one()
            } else {
                BigInt::zero()
            };
            let actual = SchubertExpr::from_class(&g, lambda.clone())
                .mul(&SchubertExpr::from_class(&g, mu.clone()))
                .unwrap()
                .integral(&g)
                .unwrap();
            assert_eq!(actual, expected, "lambda={lambda}, mu={mu}, dual={dual}");
        }
    }
}

#[test]
fn multiplication_is_commutative_associative_and_positive_in_gr_2_4() {
    // Cup product is a commutative associative ring, with Littlewood-Richardson coefficients.
    let g = Grassmannian::new(2, 4).unwrap();
    let all = partitions(&g)
        .into_iter()
        .map(|partition| SchubertExpr::from_class(&g, partition))
        .collect::<Vec<_>>();

    for a in &all {
        for b in &all {
            let ab = a.mul(b).unwrap();
            let ba = b.mul(a).unwrap();
            assert_eq!(ab, ba);
            assert!(ab.terms().values().all(|coeff| coeff >= &BigInt::zero()));

            for c in &all {
                let left = ab.mul(c).unwrap();
                let right = a.mul(&b.mul(c).unwrap()).unwrap();
                assert_eq!(left, right);
            }
        }
    }
}

#[test]
fn degree_of_gr_2_5_under_plucker_embedding_is_five() {
    // deg Gr(2,5) is integral(sigma_1^6), counted by SYT of rectangular shape (3,3).
    let g = Grassmannian::new(2, 5).unwrap();
    let hyperplane = class(&g, &[1]);

    assert_eq!(
        hyperplane.pow(6).unwrap().integral(&g).unwrap(),
        BigInt::from(5)
    );
}
