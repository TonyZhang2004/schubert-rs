//! Classical Littlewood-Richardson multiplication for Grassmannian Schubert
//! classes.
//!
//! The coefficient `c^nu_{lambda,mu}` is counted by semistandard tableaux of
//! skew shape `nu / lambda` and content `mu`.  Cells are filled in reading-word
//! order: top row to bottom row, right to left within each row.  Every prefix of
//! that word must be lattice/Yamanouchi.

use crate::{Grassmannian, Partition, Result, SchubertExpr};
use num_bigint::BigInt;
use num_traits::{One, Zero};

pub(crate) fn multiply_expr_by_class(expr: &SchubertExpr, mu: &Partition) -> Result<SchubertExpr> {
    let g = expr.grassmannian();
    let mut out = SchubertExpr::zero(g);

    if expr.is_zero() {
        return Ok(out);
    }

    for (lambda, coeff) in expr.terms() {
        for (nu, lr_coeff) in multiply_partitions(g, lambda, mu) {
            out.add_term(nu, coeff * lr_coeff);
        }
    }

    Ok(out)
}

fn multiply_partitions(
    g: &Grassmannian,
    lambda: &Partition,
    mu: &Partition,
) -> Vec<(Partition, BigInt)> {
    let target_size = lambda.size() + mu.size();
    if target_size > g.dimension() {
        return Vec::new();
    }

    let content = mu
        .parts()
        .iter()
        .map(|&part| part as usize)
        .collect::<Vec<_>>();

    partitions_containing(g, lambda, target_size)
        .into_iter()
        .filter_map(|nu| {
            let coefficient = count_tableaux(g, lambda, &nu, &content);
            (!coefficient.is_zero()).then_some((nu, coefficient))
        })
        .collect()
}

fn partitions_containing(
    g: &Grassmannian,
    lambda: &Partition,
    target_size: usize,
) -> Vec<Partition> {
    let mut builder = PartitionBuilder {
        g,
        lambda: lambda.parts(),
        target_size,
        current: Vec::with_capacity(g.k()),
        out: Vec::new(),
    };
    builder.generate(0, g.width() as u16, 0);
    builder.out
}

struct PartitionBuilder<'a> {
    g: &'a Grassmannian,
    lambda: &'a [u16],
    target_size: usize,
    current: Vec<u16>,
    out: Vec<Partition>,
}

impl PartitionBuilder<'_> {
    fn generate(&mut self, row: usize, max_part: u16, current_size: usize) {
        if row == self.g.k() {
            if current_size == self.target_size {
                self.out
                    .push(Partition::from_normalized_unchecked(self.current.clone()));
            }
            return;
        }

        let rows_left_after = self.g.k() - row - 1;
        let lower = self.lambda[row];
        for part in (lower..=max_part).rev() {
            let next_size = current_size + part as usize;
            let min_possible = next_size
                + self.lambda[row + 1..]
                    .iter()
                    .map(|&part| part as usize)
                    .sum::<usize>();
            let max_possible = next_size + rows_left_after * part as usize;

            if self.target_size < min_possible || self.target_size > max_possible {
                continue;
            }

            self.current.push(part);
            self.generate(row + 1, part, next_size);
            self.current.pop();
        }
    }
}

fn count_tableaux(
    g: &Grassmannian,
    lambda: &Partition,
    nu: &Partition,
    content: &[usize],
) -> BigInt {
    let cells = reading_word_cells(lambda, nu);
    if cells.iter().any(|&(_, col)| col >= g.width()) {
        return BigInt::zero();
    }

    if cells.len() != content.iter().sum::<usize>() {
        return BigInt::zero();
    }

    let mut state = TableauState {
        lambda: lambda.parts(),
        nu: nu.parts(),
        board: vec![vec![None; g.width()]; g.k()],
        remaining: content.to_vec(),
        prefix_counts: vec![0; content.len()],
    };

    count_fillings(&cells, 0, &mut state)
}

fn reading_word_cells(lambda: &Partition, nu: &Partition) -> Vec<(usize, usize)> {
    let mut cells = Vec::new();
    for row in 0..lambda.parts().len() {
        for col in (lambda.parts()[row] as usize..nu.parts()[row] as usize).rev() {
            cells.push((row, col));
        }
    }
    cells
}

struct TableauState<'a> {
    lambda: &'a [u16],
    nu: &'a [u16],
    board: Vec<Vec<Option<usize>>>,
    remaining: Vec<usize>,
    prefix_counts: Vec<usize>,
}

fn count_fillings(cells: &[(usize, usize)], index: usize, state: &mut TableauState<'_>) -> BigInt {
    if index == cells.len() {
        return BigInt::one();
    }

    let (row, col) = cells[index];
    let mut total = BigInt::zero();

    for value in 0..state.remaining.len() {
        if state.remaining[value] == 0 || !is_semistandard_at(state, row, col, value) {
            continue;
        }

        state.remaining[value] -= 1;
        state.prefix_counts[value] += 1;

        if is_lattice_prefix(&state.prefix_counts) {
            state.board[row][col] = Some(value);
            total += count_fillings(cells, index + 1, state);
            state.board[row][col] = None;
        }

        state.prefix_counts[value] -= 1;
        state.remaining[value] += 1;
    }

    total
}

fn is_semistandard_at(state: &TableauState<'_>, row: usize, col: usize, value: usize) -> bool {
    if col + 1 < state.nu[row] as usize && col + 1 >= state.lambda[row] as usize {
        if let Some(right) = state.board[row][col + 1] {
            if value > right {
                return false;
            }
        }
    }

    if row > 0 && col < state.nu[row - 1] as usize && col >= state.lambda[row - 1] as usize {
        if let Some(above) = state.board[row - 1][col] {
            if above >= value {
                return false;
            }
        }
    }

    true
}

fn is_lattice_prefix(prefix_counts: &[usize]) -> bool {
    prefix_counts
        .windows(2)
        .all(|counts| counts[0] >= counts[1])
}
