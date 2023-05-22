use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::Field;

pub mod mipp;
mod kzg;

// Computes `l + xr` pointwise.
fn fold_points<A: AffineRepr>(l: &[A], r: &[A], x: &A::ScalarField) -> Vec<A> {
    assert_eq!(l.len(), r.len());
    let proj: Vec<A::Group> = l.iter().zip(r).map(|(&l, &r)| (l + r * x)).collect();
    A::Group::normalize_batch(&proj)
}

// Computes `l + xr` pointwise.
fn fold_scalars<A: Field>(l: &[A], r: &[A], x: &A) -> Vec<A> {
    l.iter().zip(r).map(|(&l, &r)| (l + r * x)).collect()
}

// n = 2^m
// Folding elements V = [A1, ..., An] with scalars [x1, ..., xm] recursively m times using formula
// V = VL || VR, Vi = VL + xi * VR pointwise, V := Vi, i = 1,...,m
// results in V = Vm = [c1A1 + ... + cnAn], where ci = prod({xj | if j-th bit of i-1 is set}).
// This function computes these ci-s.
fn final_folding_exponents<F: Field>(xs: &[F]) -> Vec<F> {
    let m = xs.len();
    let mut n = 2usize.pow(m as u32);
    let mut res = vec![F::one(); n];
    for x in xs {
        n = n / 2;
        for chunk in res.rchunks_mut(n).step_by(2) {
            for elem in chunk.iter_mut() {
                *elem *= x;
            }
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_std::{test_rng, UniformRand};

    use super::*;

    #[test]
    fn test_final_folding_exponents() {
        let rng = &mut test_rng();

        let m = 10;
        let n = 2usize.pow(m);

        let xs: Vec<Fr> = (0..m).map(|_| Fr::rand(rng)).collect();
        let v: Vec<Fr> = (0..n).map(|_| Fr::rand(rng)).collect();

        let mut vi = v.clone();
        let mut mi = vi.len();
        for x in xs.iter() {
            mi /= 2;
            vi = fold_scalars(&vi[..mi],&vi[mi..], x);
        }
        let final_v = vi[0];

        let cs = final_folding_exponents(&xs);
        let multi_exp: Fr = v.iter().zip(cs.iter()).map(|(x, y)| x * y).sum();

        assert_eq!(final_v, multi_exp);
    }
}