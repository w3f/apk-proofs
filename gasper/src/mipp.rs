use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ff::{CyclotomicMultSubgroup, Field, PrimeField};
use ark_std::{test_rng, UniformRand};
use crate::{fold_msm, fold_scalars};

pub struct ProverKey<E: Pairing> {
    v: Vec<E::G2Affine>,
    w: Vec<E::G1Affine>,
}

pub struct VerifierKey<E: Pairing> {
    v: Vec<E::G2Affine>,
    w: Vec<E::G1Affine>,
}

pub struct Proof<E: Pairing> {
    a_comm_l: E::TargetField,
    a_comm_r: E::TargetField,
    b_comm_l: E::G1,
    b_comm_r: E::G1,
    c_l: E::G1,
    c_r: E::G1,
    a: E::G1,
    b: E::ScalarField,
    x: E::ScalarField,
}

fn setup<E: Pairing>(log_m: u32) -> (ProverKey<E>, VerifierKey<E>) {
    let rng = &mut test_rng();

    let m = 2usize.pow(log_m);

    let v: Vec<E::G2Affine> = (0..m).map(|_| E::G2Affine::rand(rng)).collect();
    let w: Vec<E::G1Affine> = (0..m).map(|_| E::G1Affine::rand(rng)).collect();

    let pk = ProverKey {
        v: v.clone(),
        w: w.clone(),
    };
    let vk = VerifierKey {
        v,
        w,
    };
    (pk, vk)
}

pub fn prove<E: Pairing>(pk: &ProverKey<E>, a: &[E::G1Affine], b: &[E::ScalarField]) -> Proof<E> {
    let rng = &mut test_rng();

    let m = a.len();
    assert_eq!(b.len(), m);

    let v = &pk.v;
    let w = &pk.w;

    let m1 = m / 2;

    let a_l = &a[..m1];
    let a_r = &a[m1..];
    let b_l = &b[..m1];
    let b_r = &b[m1..];

    let v_l = &v[..m1];
    let v_r = &v[m1..];
    let w_l = &w[..m1];
    let w_r = &w[m1..];

    let a_comm_l: E::TargetField = E::multi_pairing(a_r, v_l).0;
    let a_comm_r: E::TargetField = E::multi_pairing(a_l, v_r).0;
    let b_comm_l: E::G1 = VariableBaseMSM::msm(&w_r, &b_l).unwrap();
    let b_comm_r: E::G1 = VariableBaseMSM::msm(&w_l, &b_r).unwrap();
    let c_l: E::G1 = VariableBaseMSM::msm(&a_r, &b_l).unwrap();
    let c_r: E::G1 = VariableBaseMSM::msm(&a_l, &b_r).unwrap();

    let x = E::ScalarField::rand(rng);
    let x_inv = x.inverse().unwrap();

    let a1 = fold_msm(&a, &x);
    let b1 = fold_scalars(&b, &x_inv);
    // let v1 = fold_msm(&v, &x_inv);
    // let w1 = fold_msm(&w, &x);
    Proof{
        a_comm_l,
        a_comm_r,
        b_comm_l,
        b_comm_r,
        c_l,
        c_r,
        a: a1[0].into(),
        b: b1[0],
        x,
    }
}

pub fn verify<E: Pairing>(vk: &VerifierKey<E>, proof: &Proof<E>, a_comm: &E::TargetField, b_comm: &E::G1, c: &E::G1) {
    let v = &vk.v;
    let w = &vk.w;

    let a_comm_l = proof.a_comm_l;
    let a_comm_r = proof.a_comm_r;
    let b_comm_l = proof.b_comm_l;
    let b_comm_r = proof.b_comm_r;
    let c_l = proof.c_l;
    let c_r = proof.c_r;
    let a1 = proof.a;
    let b1 = proof.b;
    let x = proof.x;

    let x_inv = x.inverse().unwrap();

    let v1 = fold_msm(&v, &x_inv)[0];
    let w1 = fold_msm(&w, &x)[0];

    let a_comm_1 = a_comm_l.cyclotomic_exp(<E as Pairing>::ScalarField::into_bigint(x))
        * a_comm
        * a_comm_r.cyclotomic_exp(x_inv.into_bigint());
    assert_eq!(a_comm_1, E::pairing(&a1, &v1).0);
    let b_comm_1 = b_comm_l * x + b_comm + b_comm_r * x_inv;
    assert_eq!(b_comm_1, w1 * b1);
    let c_1 = c_l * x + c + c_r * x_inv;
    assert_eq!(c_1, a1 * b1);
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

        let log_m = 1;
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