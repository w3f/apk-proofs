use ark_ec::CurveGroup;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ff::{batch_inversion, CyclotomicMultSubgroup, Field, PrimeField};
use ark_std::{test_rng, UniformRand};

use crate::{fold_msm, fold_scalars};

pub struct ProverKey<E: Pairing> {
    log_m: u32,
    v: Vec<E::G2Affine>,
    w: Vec<E::G1Affine>,
    h1: E::G2Affine,
    h2: E::G2Affine,
}

pub struct VerifierKey<E: Pairing> {
    log_m: u32,
    v: Vec<E::G2Affine>,
    w: Vec<E::G1Affine>,
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

    xs: Vec<E::ScalarField>,
    c1: E::ScalarField,
    c2: E::ScalarField,
}

fn setup<E: Pairing>(log_m: u32) -> (ProverKey<E>, VerifierKey<E>) {
    let rng = &mut test_rng();

    let m = 2usize.pow(log_m);

    let v: Vec<E::G2Affine> = (0..m).map(|_| E::G2Affine::rand(rng)).collect();
    let w: Vec<E::G1Affine> = (0..m).map(|_| E::G1Affine::rand(rng)).collect();
    let h1 = E::G2Affine::rand(rng);
    let h2 = E::G2Affine::rand(rng);

    let pk = ProverKey {
        log_m,
        v: v.clone(),
        w: w.clone(),
        h1,
        h2,
    };
    let vk = VerifierKey {
        log_m,
        v,
        w,
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
    let h1 = (pk.h1 * c1).into_affine();
    let h2 = (pk.h2 * c2).into_affine();

    let mut m1 = m;
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    let mut v = pk.v.to_vec();
    let mut w = pk.w.to_vec();

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

    Proof {
        comm_ls,
        comm_rs,
        a: a[0].into(),
        b: b[0],
        v: v[0],
        w: w[0],
        xs,
        c1,
        c2,
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

    let mut v = vk.v.to_vec();
    let mut w = vk.w.to_vec();

    let h1 = (vk.h1 * proof.c1).into_affine();
    let h2 = (vk.h2 * proof.c2).into_affine();

    for (x, x_inv) in xs.iter().zip(xs_inv.iter()) {
        v = fold_msm(&v, &x_inv);
        w = fold_msm(&w, &x);
    }

    assert_eq!(v[0], proof.v);
    assert_eq!(w[0], proof.w);

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
    assert_eq!(comm, E::multi_pairing([a, w[0] * b, a * b], [v[0], h1, h2]).0);
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_ec::VariableBaseMSM;
    use ark_std::{test_rng, UniformRand};

    use crate::mipp::{prove, setup, verify};

    fn _test_mipp<E: Pairing>() {
        let rng = &mut test_rng();

        let log_m = 2;
        let m = 2usize.pow(log_m);

        let (pk, vk) = setup::<E>(log_m);

        // Want to prove <A, b> = b1A1 + ... + bmAm = C
        let a: Vec<E::G1Affine> = (0..m).map(|_| E::G1Affine::rand(rng)).collect();
        let b: Vec<E::ScalarField> = (0..m).map(|_| E::ScalarField::rand(rng)).collect();
        let c: E::G1 = VariableBaseMSM::msm(&a, &b).unwrap();

        let v = &pk.v;
        let w = &pk.w;

        // A_comm = <A, V> = e(A1, V1) * ... * e(Am, Vm)
        let a_comm: E::TargetField = E::multi_pairing(&a, v).0;
        // b_comm = <b, W> = b1W1 * ... + bmWm
        let b_comm: E::G1 = VariableBaseMSM::msm(w, &b).unwrap();

        let proof = prove(&pk, &a, &b);

        verify(&vk, &proof, &a_comm, &b_comm, &c);
    }

    #[test]
    fn test_mipp() {
        _test_mipp::<Bls12_381>();
    }
}