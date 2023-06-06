use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::Field;
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub mod mipp_u;
pub mod mipp_k;
pub mod sipp;

// Computes `l + xr` pointwise.
fn fold_points<A: AffineRepr>(l: &[A], r: &[A], x: &A::ScalarField) -> Vec<A> {
    assert_eq!(l.len(), r.len());
    let proj: Vec<A::Group> = ark_std::cfg_iter!(l)
        .zip(r)
        .map(|(&l, &r)| (l + r * x))
        .collect();
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

// Computes the final commitment key polynomial for the contiguous SRS (in G1).
// n = 2^m, xs = [x1, ..., xm]
// fw(X) = (1 + x_m X) (1 + x_{m-1} X^2) ... (1 + x_1 X^{2^{m-1}})
// deg(fw) = n - 1
fn compute_final_poly_for_g1<F: Field>(xs: &[F]) -> DensePolynomial<F> {
    let coeffs = final_folding_exponents(xs);
    DensePolynomial::from_coefficients_vec(coeffs)
}

fn evaluate_final_poly_for_g1<F: Field>(xs: &[F], z: &F) -> F {
    evaluate_final_poly(xs, z)
}

// Computes the final commitment key polynomial for the gapped SRS (in G2).
// n = 2^m, xs = [x1, ..., xm]
// fv(X) = (1 + x_m X^2) (1 + x_{m-1} X^4) ... (1 + x_1 X^{2^m})
// deg(fv) = 2n - 2
// fv(X) = fw(X^2)
fn compute_final_poly_for_g2<F: Field>(xs: &[F]) -> DensePolynomial<F> {
    let exps = final_folding_exponents(xs);
    // interleave with 0s
    let coeffs = exps.into_iter()
        .flat_map(|exp| [exp, F::zero()])
        .collect();
    DensePolynomial::from_coefficients_vec(coeffs)
}

fn evaluate_final_poly_for_g2<F: Field>(xs: &[F], z: &F) -> F {
    let z2 = z.square();
    evaluate_final_poly(xs, &z2)
}

// Computes (1 + x_m z)(1 + x_{m-1} z^2) ... (1 + x_1 z^{2^{m-1}}).
fn evaluate_final_poly<F: Field>(xs: &[F], z: &F) -> F {
    let mut res = F::one();
    let mut z_i = z.clone();
    for x in xs.iter().rev() {
        res *= z_i * x + F::one();
        z_i = z_i.square(); //TODO: remove extra squaring
    }
    res
}


#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};
    use crate::kzg;

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

    fn fold_key<A: AffineRepr>(a: Vec<A>, xs: &[A::ScalarField]) -> A {
        let mut a = a;
        let mut m1 = a.len();
        for x in xs.iter() {
            m1 /= 2;
            a = fold_points(&a[..m1], &a[m1..], x);
        }
        assert_eq!(a.len(), 1);
        a[0]
    }

    #[test]
    fn test_final_ck_polynomial_for_g1() {
        let rng = &mut test_rng();

        let log_m = 8;
        let m = 2usize.pow(log_m);
        let xs: Vec<Fr> = (0..log_m).map(|_| Fr::rand(rng)).collect();

        let fw = compute_final_poly_for_g1(&xs);
        assert_eq!(fw.degree(), m - 1);

        let kzg_ck = {
            let g = G1Projective::rand(rng);
            let h = G2Projective::rand(rng);
            let (kzg_ck, _) = kzg::setup_g1::<Bls12_381>(m, g, h);
            kzg_ck
        };
        let w = kzg_ck.clone();
        let final_w = fold_key(w, &xs);
        let fw_comm = kzg::commit_g1::<Bls12_381>(&kzg_ck, &fw);
        assert_eq!(fw_comm, final_w);

        let z = Fr::rand(rng);
        let fw_at_z_1 = fw.evaluate(&z);
        let fw_at_z_2 = evaluate_final_poly_for_g1(&xs, &z);
        assert_eq!(fw_at_z_1, fw_at_z_2);
    }

    #[test]
    fn test_final_ck_polynomial_for_g2() {
        let rng = &mut test_rng();

        let log_m = 8;
        let m = 2usize.pow(log_m);
        let xs: Vec<Fr> = (0..log_m).map(|_| Fr::rand(rng)).collect();

        let fv = compute_final_poly_for_g2(&xs);
        assert_eq!(fv.degree(), 2 * m - 2);

        let kzg_ck = {
            let g = G1Projective::rand(rng);
            let h = G2Projective::rand(rng);
            let (kzg_ck, _) = kzg::setup_g2::<Bls12_381>(2 * m - 1, g, h);
            kzg_ck
        };
        let v = kzg_ck.iter().cloned().step_by(2).collect();
        let final_v = fold_key(v, &xs);
        let fv_comm = kzg::commit_g2::<Bls12_381>(&kzg_ck, &fv);
        assert_eq!(fv_comm, final_v);

        let z = Fr::rand(rng);
        let fv_at_z_1 = fv.evaluate(&z);
        let fv_at_z_2 = evaluate_final_poly_for_g2(&xs, &z);
        assert_eq!(fv_at_z_1, fv_at_z_2);
    }
}