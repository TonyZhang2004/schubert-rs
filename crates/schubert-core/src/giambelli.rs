use crate::{Partition, Result, SchubertExpr};
use num_bigint::BigInt;
use num_traits::One;

pub(crate) fn multiply_expr_by_class(expr: &SchubertExpr, mu: &Partition) -> Result<SchubertExpr> {
    let g = expr.grassmannian();
    let mut out = SchubertExpr::zero(g);

    for permutation in permutations(g.k()) {
        let mut specials = Vec::with_capacity(g.k());
        let mut is_zero_term = false;

        for (i, &p_i) in permutation.iter().enumerate() {
            let index = mu.parts()[i] as isize + p_i as isize - i as isize;
            if index < 0 || index as usize > g.width() {
                is_zero_term = true;
                break;
            }
            specials.push(index);
        }

        if is_zero_term {
            continue;
        }

        let mut term = expr.clone();
        for special in specials {
            term = term.multiply_by_special(special)?;
        }

        if permutation_sign(&permutation) < 0 {
            term = term.scalar_mul(&BigInt::from(-1));
        } else {
            term = term.scalar_mul(&BigInt::one());
        }

        out = out.checked_add(&term)?;
    }

    Ok(out)
}

fn permutations(n: usize) -> Vec<Vec<usize>> {
    fn backtrack(start: usize, values: &mut [usize], out: &mut Vec<Vec<usize>>) {
        if start == values.len() {
            out.push(values.to_vec());
            return;
        }

        for index in start..values.len() {
            values.swap(start, index);
            backtrack(start + 1, values, out);
            values.swap(start, index);
        }
    }

    let mut values = (0..n).collect::<Vec<_>>();
    let mut out = Vec::new();
    backtrack(0, &mut values, &mut out);
    out
}

fn permutation_sign(permutation: &[usize]) -> i8 {
    let inversions = permutation
        .iter()
        .enumerate()
        .map(|(i, &left)| {
            permutation[i + 1..]
                .iter()
                .filter(|&&right| left > right)
                .count()
        })
        .sum::<usize>();

    if inversions % 2 == 0 {
        1
    } else {
        -1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn permutation_sign_uses_inversion_parity() {
        assert_eq!(permutation_sign(&[0, 1, 2]), 1);
        assert_eq!(permutation_sign(&[1, 0, 2]), -1);
        assert_eq!(permutation_sign(&[2, 1, 0]), -1);
    }
}
