pub mod mipp;

use ark_ec::AffineRepr;
use ark_ff::Field;

// TODO: MSM
fn fold_msm<A: AffineRepr>(a: &[A], x: &A::ScalarField) -> Vec<A> {
    let m = a.len() / 2;
    let (l, r) = a.split_at(m);
    l.iter().zip(r).map(|(&l, &r)| (l + r * x).into()).collect()
}

fn fold_scalars<A: Field>(a: &[A], x: &A) -> Vec<A> {
    let m = a.len() / 2;
    let (l, r) = a.split_at(m);
    l.iter().zip(r).map(|(&l, &r)| (l + r * x)).collect()
}