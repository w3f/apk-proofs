use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ff::{batch_inversion, CyclotomicMultSubgroup, Field, PrimeField};
use ark_std::{test_rng, UniformRand};
use crate::{fold_msm, fold_scalars};

pub struct ProverKey<E: Pairing> {
    log_m: u32,
    v: Vec<E::G2Affine>,
    w: Vec<E::G1Affine>,
}

pub struct VerifierKey<E: Pairing> {
    log_m: u32,
    v: Vec<E::G2Affine>,
    w: Vec<E::G1Affine>,
}

pub struct Proof<E: Pairing> {
    a_comm_ls: Vec<E::TargetField>,
    a_comm_rs: Vec<E::TargetField>,
    b_comm_ls: Vec<E::G1>,
    b_comm_rs: Vec<E::G1>,
    c_ls: Vec<E::G1>,
    c_rs: Vec<E::G1>,
    a: E::G1,
    b: E::ScalarField,
    xs: Vec<E::ScalarField>,
}

fn setup<E: Pairing>(log_m: u32) -> (ProverKey<E>, VerifierKey<E>) {
    let rng = &mut test_rng();

    let m = 2usize.pow(log_m);

    let v: Vec<E::G2Affine> = (0..m).map(|_| E::G2Affine::rand(rng)).collect();
    let w: Vec<E::G1Affine> = (0..m).map(|_| E::G1Affine::rand(rng)).collect();

    let pk = ProverKey {
        log_m,
        v: v.clone(),
        w: w.clone(),
    };
    let vk = VerifierKey {
        log_m,
        v,
        w,
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

    let mut m1 = m;
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    let mut v = pk.v.to_vec();
    let mut w = pk.w.to_vec();

    let mut a_comm_ls: Vec<E::TargetField> = Vec::with_capacity(log_m);
    let mut a_comm_rs: Vec<E::TargetField> = Vec::with_capacity(log_m);
    let mut b_comm_ls: Vec<E::G1> = Vec::with_capacity(log_m);
    let mut b_comm_rs: Vec<E::G1> = Vec::with_capacity(log_m);
    let mut c_ls: Vec<E::G1> = Vec::with_capacity(log_m);
    let mut c_rs: Vec<E::G1> = Vec::with_capacity(log_m);

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

        let a_comm_l = E::multi_pairing(a_r, v_l).0;
        let a_comm_r = E::multi_pairing(a_l, v_r).0;
        let b_comm_l = VariableBaseMSM::msm(&w_r, &b_l).unwrap();
        let b_comm_r = VariableBaseMSM::msm(&w_l, &b_r).unwrap();
        let c_l = VariableBaseMSM::msm(&a_r, &b_l).unwrap();
        let c_r = VariableBaseMSM::msm(&a_l, &b_r).unwrap();

        a_comm_ls.push(a_comm_l);
        a_comm_rs.push(a_comm_r);
        b_comm_ls.push(b_comm_l);
        b_comm_rs.push(b_comm_r);
        c_ls.push(c_l);
        c_rs.push(c_r);

        let x_inv = x.inverse().unwrap();

        a = fold_msm(&a, &x);
        b = fold_scalars(&b, &x_inv);
        v = fold_msm(&v, &x_inv);
        w = fold_msm(&w, &x);
    }

    Proof {
        a_comm_ls,
        a_comm_rs,
        b_comm_ls,
        b_comm_rs,
        c_ls,
        c_rs,
        a: a[0].into(),
        b: b[0],
        xs,
    }
}

pub fn verify<E: Pairing>(vk: &VerifierKey<E>, proof: &Proof<E>, a_comm: &E::TargetField, b_comm: &E::G1, c: &E::G1) {
    let a = proof.a;
    let b = proof.b;
    let a_comm_ls = proof.a_comm_ls.to_vec();
    let a_comm_rs = proof.a_comm_rs.to_vec();
    let b_comm_ls = proof.b_comm_ls.to_vec();
    let b_comm_rs = proof.b_comm_rs.to_vec();
    let c_ls = proof.c_ls.to_vec();
    let c_rs = proof.c_rs.to_vec();

    let xs = proof.xs.to_vec();
    let mut xs_inv = xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let mut v = vk.v.to_vec();
    let mut w = vk.w.to_vec();

    for (x, x_inv) in xs.iter().zip(xs_inv.iter()) {
        v = fold_msm(&v, &x_inv);
        w = fold_msm(&w, &x);
    }

    let mut a_comm = a_comm.clone();
    for (((x, x_inv), a_comm_l), a_comm_r) in xs.iter()
        .zip(xs_inv.iter())
        .zip(a_comm_ls)
        .zip(a_comm_rs) {
        a_comm = a_comm_l.cyclotomic_exp(x.into_bigint())
            * a_comm
            * a_comm_r.cyclotomic_exp(x_inv.into_bigint());
    }
    assert_eq!(a_comm, E::pairing(&a, &v[0]).0);

    let mut b_comm = b_comm.clone();
    for (((x, x_inv), b_comm_l), b_comm_r) in xs.iter()
        .zip(xs_inv.iter())
        .zip(b_comm_ls)
        .zip(b_comm_rs) {
        b_comm = b_comm_l * x + b_comm + b_comm_r * x_inv;
    }
    assert_eq!(b_comm, w[0] * b);

    let mut c = c.clone();
    for (((x, x_inv), c_l),c_r) in xs.iter()
        .zip(xs_inv.iter())
        .zip(c_ls)
        .zip(c_rs) {
        c = c_l * x + c + c_r * x_inv;
    }

    assert_eq!(c, a * b);
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

        let log_m = 10;
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