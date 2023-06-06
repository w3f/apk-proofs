use ark_ff::Field;

pub mod ipa;
mod kzg;

fn powers<F: Field>(x: F, n: usize) -> Vec<F> {
    let x0 = F::one();
    ark_std::iter::successors(Some(x0), move |xi| Some(x * xi))
        .take(n)
        .collect()
}