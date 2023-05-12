use ark_ec::AffineRepr;
use ark_ff::{FftField, Field};
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};

pub mod mipp;
mod kzg;

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

// TODO: no need in ffts for that
fn final_ck_polynomial<F: FftField>(xs: &[F]) -> DensePolynomial<F> {
    let mut f = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    for (i, &x) in xs.into_iter().rev().enumerate() {
        let n = 2usize.pow((i + 1) as u32);
        let fi = SparsePolynomial::from_coefficients_vec(vec![(0, F::one()), (n, x)]);
        let fi: DensePolynomial<F> = fi.into();
        f = &f * &fi;
    }
    f
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective};
    use ark_ec::{CurveGroup, VariableBaseMSM};
    use ark_poly::Polynomial;
    use ark_std::test_rng;
    use ark_std::UniformRand;
    use fflonk::pcs::kzg::urs::URS;

    use crate::{final_ck_polynomial, fold_msm};

    #[test]
    fn final_commitment_keys() {
        let rng = &mut test_rng();

        let log_m = 10;
        let m = 2usize.pow(log_m);
        let gts= URS::<Bls12_381>::generate(2 * m - 1, 0, rng).powers_in_g1;

        let ck: Vec<G1Affine> = gts.iter().cloned().step_by(2).collect();
        let xs: Vec<Fr> = (0..log_m).map(|_| Fr::rand(rng)).collect();

        let mut ck_folded = ck;
        for x in xs.iter() {
            ck_folded = fold_msm(&ck_folded, x);
        }
        assert_eq!(ck_folded.len(), 1);
        let final_ck = ck_folded[0];

        let final_ck_poly = final_ck_polynomial(&xs);
        assert_eq!(final_ck_poly.degree(), 2 * m - 2);

        let final_ck_poly_comm: G1Projective = VariableBaseMSM::msm(&gts, &final_ck_poly.coeffs).unwrap();
        let final_ck_poly_comm = final_ck_poly_comm.into_affine();

        assert_eq!(final_ck, final_ck_poly_comm);
    }
}