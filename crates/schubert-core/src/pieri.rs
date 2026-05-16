use crate::{Grassmannian, Partition, Result};

pub(crate) fn multiply_partition_by_special(
    g: &Grassmannian,
    lambda: &Partition,
    r: usize,
) -> Result<Vec<Partition>> {
    lambda.validate_for(g)?;

    if r == 0 {
        return Ok(vec![lambda.clone()]);
    }

    if r > g.width() {
        return Ok(Vec::new());
    }

    let mut out = Vec::new();
    let mut current = vec![0; g.k()];
    generate_horizontal_strips(g, lambda.parts(), r, 0, &mut current, &mut out)?;
    Ok(out)
}

fn generate_horizontal_strips(
    g: &Grassmannian,
    lambda: &[u16],
    remaining: usize,
    row: usize,
    current: &mut [u16],
    out: &mut Vec<Partition>,
) -> Result<()> {
    if row == g.k() {
        if remaining == 0 {
            out.push(Partition::from_normalized_unchecked(current.to_vec()));
        }
        return Ok(());
    }

    let lower = lambda[row] as usize;
    let upper = if row == 0 {
        g.width()
    } else {
        g.width().min(lambda[row - 1] as usize)
    };

    for nu_i in lower..=upper {
        let added = nu_i - lower;
        if added > remaining {
            break;
        }

        current[row] = nu_i as u16;
        generate_horizontal_strips(g, lambda, remaining - added, row + 1, current, out)?;
    }

    Ok(())
}
