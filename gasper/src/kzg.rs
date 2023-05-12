use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_ff::{One, Zero};
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial;
use ark_std::{test_rng, UniformRand};
use fflonk::utils::ec::single_base_msm;

pub(crate) struct VerifierKeyG1<E: Pairing> {
    g1: E::G1Affine,
    g2: E::G2Prepared,
    tg2: E::G2Prepared,
}

pub(crate) struct VerifierKeyG2<E: Pairing> {
    g1: E::G1Affine,
    g2: E::G2Affine,
    tg1: E::G1,
}

pub(crate) fn setup_g1<E: Pairing>(n: usize, g1: E::G1, g2: E::G2) -> (Vec<E::G1Affine>, VerifierKeyG1<E>) {
    let rng = &mut test_rng();
    let t = E::ScalarField::rand(rng);
    let ts: Vec<E::ScalarField> = {
        let t0 = E::ScalarField::one();
        ark_std::iter::successors(Some(t0), move |ti| Some(t * ti)).take(n).collect()
    };
    let tsg1 = single_base_msm(&ts, g1);
    let vk = VerifierKeyG1 {
        g1: g1.into(),
        g2: g2.into_affine().into(),
        tg2: (g2 * t).into_affine().into(),
    };
    (tsg1, vk)
}

pub(crate) fn setup_g2<E: Pairing>(n: usize, g1: E::G1, g2: E::G2) -> (Vec<E::G2Affine>, VerifierKeyG2<E>) {
    let rng = &mut test_rng();
    let t = E::ScalarField::rand(rng);
    let ts: Vec<E::ScalarField> = {
        let t0 = E::ScalarField::one();
        ark_std::iter::successors(Some(t0), move |ti| Some(t * ti)).take(n).collect()
    };
    let tsg2 = single_base_msm(&ts, g2);
    let vk = VerifierKeyG2 {
        g1: g1.into(),
        g2: g2.into(),
        tg1: g1 * t,
    };
    (tsg2, vk)
}

pub(crate) fn commit_g1<E: Pairing>(ck: &[E::G1Affine], p: &DensePolynomial<E::ScalarField>) -> E::G1Affine {
    assert!(p.coeffs.len() <= ck.len());
    let c: E::G1 = VariableBaseMSM::msm_unchecked(ck, &p.coeffs);
    c.into_affine()
}

pub(crate) fn commit_g2<E: Pairing>(ck: &[E::G2Affine], p: &DensePolynomial<E::ScalarField>) -> E::G2Affine {
    assert!(p.coeffs.len() <= ck.len());
    let c: E::G2 = VariableBaseMSM::msm_unchecked(ck, &p.coeffs);
    c.into_affine()
}

pub(crate) fn open_g1<E: Pairing>(
    ck: &[E::G1Affine],
    p: &DensePolynomial<E::ScalarField>,
    x: E::ScalarField,
) -> E::G1Affine {
    let v = DensePolynomial::from_coefficients_vec(vec![-x, E::ScalarField::one()]);
    let q = p / &v;
    commit_g1::<E>(&ck, &q)
}

pub(crate) fn open_g2<E: Pairing>(
    ck: &[E::G2Affine],
    p: &DensePolynomial<E::ScalarField>,
    x: E::ScalarField,
) -> E::G2Affine {
    let v = DensePolynomial::from_coefficients_vec(vec![-x, E::ScalarField::one()]);
    let q = p / &v;
    commit_g2::<E>(&ck, &q)
}

pub(crate) fn verify_g1<E: Pairing>(
    vk: &VerifierKeyG1<E>,
    c: E::G1Affine,
    x: E::ScalarField,
    y: E::ScalarField,
    proof: E::G1Affine,
) -> bool {
    let acc = vk.g1 * y - (proof * x + c);
    let acc = acc.into_affine();
    E::multi_pairing([acc, proof], [vk.g2.clone(), vk.tg2.clone()]).is_zero()
}

// TODO: bench
// I think this formula is beneficial for the G2 case.
pub(crate) fn verify_g2<E: Pairing>(
    vk: &VerifierKeyG2<E>,
    c: E::G2Affine,
    x: E::ScalarField,
    y: E::ScalarField,
    proof: E::G2Affine,
) -> bool {
    let g2_val = (vk.g2 * y - c).into_affine();
    let g1_val = (vk.tg1 - vk.g1 * x).into_affine();
    E::multi_pairing([vk.g1, g1_val], [g2_val, proof]).is_zero()
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Affine, G2Projective};
    use ark_poly::DenseUVPolynomial;
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};

    use super::*;

    #[test]
    fn test_kzg_g1() {
        let rng = &mut test_rng();

        let n = 16;

        let g1 = G1Projective::rand(rng);
        let g2 = G2Projective::rand(rng);

        let (ck, vk) = setup_g1::<Bls12_381>(n, g1, g2);

        let p = DensePolynomial::<Fr>::rand(n - 1, rng);
        let x = Fr::rand(rng);
        let y = p.evaluate(&x);

        let c: G1Affine = commit_g1::<Bls12_381>(&ck, &p);

        let proof = open_g1::<Bls12_381>(&ck, &p, x);

        assert!(verify_g1(&vk, c, x, y, proof));
    }

    #[test]
    fn test_kzg_g2() {
        let rng = &mut test_rng();

        let n = 16;

        let g1 = G1Projective::rand(rng);
        let g2 = G2Projective::rand(rng);

        let (ck, vk) = setup_g2::<Bls12_381>(n, g1, g2);

        let p = DensePolynomial::<Fr>::rand(n - 1, rng);
        let x = Fr::rand(rng);
        let y = p.evaluate(&x);

        let c: G2Affine = commit_g2::<Bls12_381>(&ck, &p);

        let proof = open_g2::<Bls12_381>(&ck, &p, x);

        assert!(verify_g2(&vk, c, x, y, proof));
    }
}