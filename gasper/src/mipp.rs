use ark_ec::CurveGroup;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ff::{batch_inversion, CyclotomicMultSubgroup, Field, PrimeField};
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_std::{test_rng, UniformRand};

use crate::{final_folding_exponents, fold_msm, fold_scalars, kzg};

pub struct ProverKey<E: Pairing> {
    log_m: u32,
    ck_g1: Vec<E::G1Affine>,
    ck_g2: Vec<E::G2Affine>,
    h1: E::G2Affine,
    h2: E::G2Affine,
}

pub struct VerifierKey<E: Pairing> {
    vk_g1: kzg::VerifierKeyG1<E>,
    vk_g2: kzg::VerifierKeyG2<E>,
    h1: E::G2Affine,
    h2: E::G2Affine,
}

pub struct Proof<E: Pairing> {
    comm_ls: Vec<E::TargetField>,
    comm_rs: Vec<E::TargetField>,
    a: E::G1,
    b: E::ScalarField,
    v: E::G2Affine,
    w: E::G1Affine,
    kzg_g1: E::G1Affine,
    kzg_g2: E::G2Affine,

    xs: Vec<E::ScalarField>,
    c1: E::ScalarField,
    c2: E::ScalarField,
    z: E::ScalarField,
}

pub fn setup<E: Pairing>(log_m: u32) -> (ProverKey<E>, VerifierKey<E>) {
    let rng = &mut test_rng();

    let m = 2usize.pow(log_m);
    let g = E::G1::rand(rng);
    let h = E::G2::rand(rng);

    let (ck_g1, vk_g1) = kzg::setup_g1::<E>(2 * m - 1, g, h);
    let (ck_g2, vk_g2) = kzg::setup_g2::<E>(2 * m - 1, g, h);

    let h1 = E::G2Affine::rand(rng);
    let h2 = E::G2Affine::rand(rng);

    let pk = ProverKey {
        log_m,
        ck_g1,
        ck_g2,
        h1,
        h2,
    };
    let vk = VerifierKey {
        vk_g1,
        vk_g2,
        h1,
        h2,
    };
    (pk, vk)
}

pub fn prove<E: Pairing>(pk: &ProverKey<E>, a: &[E::G1Affine], b: &[E::ScalarField]) -> Proof<E> {
    let rng = &mut test_rng();

    let log_m = pk.log_m as usize;
    let m = 2usize.pow(log_m as u32);
    assert_eq!(a.len(), m);
    assert_eq!(b.len(), m);

    let xs: Vec<E::ScalarField> = (0..log_m).map(|_| E::ScalarField::rand(rng)).collect();
    let c1 = E::ScalarField::rand(rng);
    let c2 = E::ScalarField::rand(rng);
    let z = E::ScalarField::rand(rng);

    let h1 = (pk.h1 * c1).into_affine();
    let h2 = (pk.h2 * c2).into_affine();

    let mut m1 = m;
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    let mut v: Vec<E::G2Affine> = pk.ck_g2.iter().cloned().step_by(2).collect();
    let mut w: Vec<E::G1Affine> = pk.ck_g1.iter().cloned().step_by(2).collect();

    let mut comm_ls: Vec<E::TargetField> = Vec::with_capacity(log_m);
    let mut comm_rs: Vec<E::TargetField> = Vec::with_capacity(log_m);

    for x in xs.iter() {
        m1 /= 2;

        let a_l = &a[..m1];
        let a_r = &a[m1..];
        let b_l = &b[..m1];
        let b_r = &b[m1..];

        let v_l = &v[..m1];
        let v_r = &v[m1..];
        let w_l = &w[..m1];
        let w_r = &w[m1..];

        let b_comm_l: E::G1 = VariableBaseMSM::msm(&w_r, &b_l).unwrap();
        let b_comm_r: E::G1 = VariableBaseMSM::msm(&w_l, &b_r).unwrap();
        let c_l: E::G1 = VariableBaseMSM::msm(&a_r, &b_l).unwrap();
        let c_r: E::G1 = VariableBaseMSM::msm(&a_l, &b_r).unwrap();

        let b_comm_l = b_comm_l.into_affine();
        let b_comm_r = b_comm_r.into_affine();
        let c_l = c_l.into_affine();
        let c_r = c_r.into_affine();

        let vals_l = [a_r, &[b_comm_l, c_l]].concat();
        let vals_r = [a_l, &[b_comm_r, c_r]].concat();
        let keys_l = [v_l, &[h1, h2]].concat();
        let keys_r = [v_r, &[h1, h2]].concat();

        let comm_l = E::multi_pairing(vals_l, keys_l).0;
        let comm_r = E::multi_pairing(vals_r, keys_r).0;

        comm_ls.push(comm_l);
        comm_rs.push(comm_r);

        let x_inv = x.inverse().unwrap();

        a = fold_msm(&a, &x);
        b = fold_scalars(&b, &x_inv);
        v = fold_msm(&v, &x_inv);
        w = fold_msm(&w, &x);
    }

    assert_eq!(a.len(), 1);
    assert_eq!(b.len(), 1);
    assert_eq!(v.len(), 1);
    assert_eq!(w.len(), 1);

    let mut xs_inv = xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let f_v = final_ck_polynomial(&xs_inv);
    let f_w = final_ck_polynomial(&xs);
    assert_eq!(kzg::commit_g2::<E>(&pk.ck_g2, &f_v), v[0]);
    assert_eq!(kzg::commit_g1::<E>(&pk.ck_g1, &f_w), w[0]);
    let kzg_g2 = kzg::open_g2::<E>(&pk.ck_g2, &f_v, z);
    let kzg_g1 = kzg::open_g1::<E>(&pk.ck_g1, &f_w, z);

    Proof {
        comm_ls,
        comm_rs,
        a: a[0].into(),
        b: b[0],
        v: v[0],
        w: w[0],
        kzg_g1,
        kzg_g2,
        xs,
        c1,
        c2,
        z,
    }
}

pub fn verify<E: Pairing>(vk: &VerifierKey<E>, proof: &Proof<E>, a_comm: &E::TargetField, b_comm: &E::G1, c: &E::G1) {
    let a = proof.a;
    let b = proof.b;
    let comm_ls = proof.comm_ls.to_vec();
    let comm_rs = proof.comm_rs.to_vec();

    let xs = proof.xs.to_vec();
    let mut xs_inv = xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let h1 = (vk.h1 * proof.c1).into_affine();
    let h2 = (vk.h2 * proof.c2).into_affine();

    let z = proof.z;
    let v = proof.v;
    let w = proof.w;

    let f_v_z = final_ck_polynomial(&xs_inv).evaluate(&z);
    let f_w_z = final_ck_polynomial(&xs).evaluate(&z);

    assert!(kzg::verify_g2(&vk.vk_g2, v, proof.z, f_v_z, proof.kzg_g2));
    assert!(kzg::verify_g1(&vk.vk_g1, w, proof.z, f_w_z, proof.kzg_g1));

    let extra = E::multi_pairing([b_comm, c], [h1, h2]).0;

    let mut comm = extra * a_comm;
    for (((x, x_inv), a_comm_l), a_comm_r) in xs.iter()
        .zip(xs_inv.iter())
        .zip(comm_ls)
        .zip(comm_rs) {
        comm = a_comm_l.cyclotomic_exp(x.into_bigint())
            * comm
            * a_comm_r.cyclotomic_exp(x_inv.into_bigint());
    }
    assert_eq!(comm, E::multi_pairing([a, w * b, a * b], [v, h1, h2]).0);
}

fn final_ck_polynomial<F: Field>(xs: &[F]) -> DensePolynomial<F> {
    let exps = final_folding_exponents(xs);
    // interleave with 0s
    let coeffs = exps.into_iter()
        .flat_map(|exp| [exp, F::zero()])
        .collect();
    DensePolynomial::from_coefficients_vec(coeffs)
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective};
    use ark_poly::univariate::SparsePolynomial;
    use ark_std::{test_rng, UniformRand};
    use fflonk::pcs::kzg::urs::URS;

    use super::*;

    fn _test_mipp<E: Pairing>() {
        let rng = &mut test_rng();

        let log_m = 2;
        let m = 2usize.pow(log_m);

        let (pk, vk) = setup::<E>(log_m);

        // Want to prove <A, b> = b1A1 + ... + bmAm = C
        let a: Vec<E::G1Affine> = (0..m).map(|_| E::G1Affine::rand(rng)).collect();
        let b: Vec<E::ScalarField> = (0..m).map(|_| E::ScalarField::rand(rng)).collect();
        let c: E::G1 = VariableBaseMSM::msm(&a, &b).unwrap();

        let v: Vec<E::G2Affine> = pk.ck_g2.iter().cloned().step_by(2).collect();
        let w: Vec<E::G1Affine> = pk.ck_g1.iter().cloned().step_by(2).collect();

        // A_comm = <A, V> = e(A1, V1) * ... * e(Am, Vm)
        let a_comm: E::TargetField = E::multi_pairing(&a, v).0;
        // b_comm = <b, W> = b1W1 * ... + bmWm
        let b_comm: E::G1 = VariableBaseMSM::msm(&w, &b).unwrap();

        let proof = prove(&pk, &a, &b);

        verify(&vk, &proof, &a_comm, &b_comm, &c);
    }

    #[test]
    fn test_mipp() {
        _test_mipp::<Bls12_381>();
    }

    // Computes final commitment key polynomial using the formula from
    // https://eprint.iacr.org/2019/1177.pdf, section 5.2
    fn naive_final_ck_polynomial<F: Field>(xs: &[F]) -> DensePolynomial<F> {
        let mut f = DensePolynomial::from_coefficients_vec(vec![F::one()]);
        for (i, &x) in xs.into_iter().rev().enumerate() {
            let n = 2usize.pow((i + 1) as u32);
            let fi = SparsePolynomial::from_coefficients_vec(vec![(0, F::one()), (n, x)]);
            f = f.naive_mul(&fi.into());
        }
        f
    }

    #[test]
    fn test_final_ck_polynomial() {
        let rng = &mut test_rng();

        let log_m = 2;
        let m = 2usize.pow(log_m);

        let xs: Vec<Fr> = (0..m).map(|_| Fr::rand(rng)).collect();

        let res = final_ck_polynomial(&xs);
        let expected = naive_final_ck_polynomial(&xs);
        assert_eq!(res, expected);
    }

    #[test]
    fn final_commitment_keys() {
        let rng = &mut test_rng();

        let log_m = 10;
        let m = 2usize.pow(log_m);
        let gts = URS::<Bls12_381>::generate(2 * m - 1, 0, rng).powers_in_g1;

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