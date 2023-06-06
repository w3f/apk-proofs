use ark_ec::CurveGroup;
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::VariableBaseMSM;
use ark_ff::{batch_inversion, Field};
use ark_std::{end_timer, start_timer, test_rng, UniformRand};
use ark_std::rand::Rng;

use crate::ipa::{compute_final_poly_for_g2, evaluate_final_poly_for_g1, evaluate_final_poly_for_g2, final_folding_exponents, fold_points, fold_scalars};
use crate::{kzg, powers};


// Inner product argument for <A, pow(b)> = A1 + bA2 + ... + b^(n-1)An = C, given the commitment
// A_comm = [A, V] = e(A1, V1) * ... * e(An, Vn), where Vi = s^(2i-2)H in G2, and
// a public scalar b

// { (A_comm, b, C; A) | C = <A, pow(b)>, A_comm = [A, V] }


pub struct ProverKey<E: Pairing> {
    log_n: u32,
    // 2n-1 points in G2: H, sH, ..., s^(2n-2)H
    pub ck_g2: Vec<E::G2Affine>,
    h1: E::G2Affine,
}

pub struct VerifierKey<E: Pairing> {
    kzg_vk_g2: kzg::VerifierKeyG2<E>,
    h1: E::G2Affine,
}

pub struct Proof<E: Pairing> {
    // Left cross-commitments
    l_comms: Vec<PairingOutput<E>>,
    // Right cross-commitments
    r_comms: Vec<PairingOutput<E>>,
    // Final folded point
    a_final: E::G1Affine,
    // Final commitment key in G2
    v_final: E::G2Affine,
    // Final commitment key argument proof for the key in G2
    kzg_proof_g2: E::G2Affine,

    challenges: Challenges<E>,
}

struct Challenges<E: Pairing> {
    xs: Vec<E::ScalarField>,
    c1: E::ScalarField,
    z: E::ScalarField,
}

impl<E: Pairing> Challenges<E> {
    fn new(log_n: usize) -> Self {
        let rng = &mut test_rng();
        Challenges {
            xs: (0..log_n).map(|_| E::ScalarField::from(rng.gen::<u128>())).collect(),
            c1: E::ScalarField::rand(rng),
            z: E::ScalarField::rand(rng),
        }
    }
}

pub fn setup<E: Pairing>(log_n: u32) -> (ProverKey<E>, VerifierKey<E>) {
    let rng = &mut test_rng();

    let n = 2usize.pow(log_n);
    let g = E::G1::rand(rng);
    let h = E::G2::rand(rng);

    let (ck_g2, kzg_vk_g2) = kzg::setup_g2::<E>(2 * n - 1, g, h);

    let h1 = E::G2Affine::rand(rng);

    let pk = ProverKey {
        log_n,
        ck_g2,
        h1,
    };

    let vk = VerifierKey {
        kzg_vk_g2,
        h1,
    };

    (pk, vk)
}

pub fn prove_for_powers<E: Pairing>(pk: &ProverKey<E>, a: &[E::G1Affine], b: E::ScalarField) -> Proof<E> {
    let powers = powers(b, a.len());
    prove_unstructured(pk, a, &powers)
}

pub fn prove_unstructured<E: Pairing>(pk: &ProverKey<E>, a: &[E::G1Affine], b: &[E::ScalarField]) -> Proof<E> {
    let log_n = pk.log_n as usize;
    let n = 2usize.pow(log_n as u32);

    // TODO: Fiat-Shamir
    let challenges = Challenges::<E>::new(log_n);

    let h1 = (pk.h1 * challenges.c1).into_affine();

    let mut n1 = n;
    let mut a_folded = a.to_vec();
    let mut b_folded = b.to_vec();
    let mut v_folded: Vec<E::G2Affine> = pk.ck_g2.iter().cloned().step_by(2).collect();
    assert_eq!(a_folded.len(), n);
    assert_eq!(b_folded.len(), n);
    assert_eq!(v_folded.len(), n);

    let mut l_comms = Vec::<PairingOutput<E>>::with_capacity(log_n);
    let mut r_comms = Vec::<PairingOutput<E>>::with_capacity(log_n);

    let t_gipa = start_timer!(|| "GIPA");
    for x in challenges.xs.iter() {
        let t_round = start_timer!(|| format!("ROUND: m = {}", n1));
        n1 /= 2;

        let al = &a_folded[..n1];
        let ar = &a_folded[n1..];
        let bl = &b_folded[..n1];
        let br = &b_folded[n1..];
        let vl = &v_folded[..n1];
        let vr = &v_folded[n1..];

        let t_msm = start_timer!(|| format!("2 x {}-msm in G1", n1));
        let cl = E::G1::msm(ar, bl).unwrap();
        let cr = E::G1::msm(al, br).unwrap();
        end_timer!(t_msm);

        // TODO: batch conversion to affine
        let cl = cl.into_affine();
        let cr = cr.into_affine();

        let l_vals = [ar, &[cl]].concat();
        let r_vals = [al, &[cr]].concat();
        let l_leys = [vl, &[h1]].concat();
        let r_keys = [vr, &[h1]].concat();

        let t_multipairing = start_timer!(|| format!("2 x {}-multipairing", n1 + 2));
        let l_comm = E::multi_pairing(l_vals, l_leys);
        let r_comm = E::multi_pairing(r_vals, r_keys);
        end_timer!(t_multipairing);
        l_comms.push(l_comm);
        r_comms.push(r_comm);

        let x_inv = x.inverse().unwrap();

        let t_folding = start_timer!(|| format!("2 x {}-folding in G1 + {}-folding in G2", n1, n1));
        a_folded = fold_points(al, ar, &x_inv);
        b_folded = fold_scalars(bl, br, &x);
        v_folded = fold_points(vl, vr, &x);
        end_timer!(t_folding);

        end_timer!(t_round);
    }
    end_timer!(t_gipa);

    assert_eq!(a_folded.len(), 1);
    assert_eq!(v_folded.len(), 1);
    let a_final = a_folded[0].into();
    let v_final = v_folded[0];

    let mut xs_inv = challenges.xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let f_v = compute_final_poly_for_g2(&challenges.xs);
    let t_kzg = start_timer!(|| format!("{}-msm in G2", 2 * n));
    let kzg_proof_g2 = kzg::open_g2::<E>(&pk.ck_g2, &f_v, challenges.z);
    end_timer!(t_kzg);

    Proof {
        l_comms,
        r_comms,
        a_final,
        v_final,
        kzg_proof_g2,
        challenges,
    }
}

pub fn verify_for_powers<E: Pairing>(vk: &VerifierKey<E>, proof: &Proof<E>, a_comm: &PairingOutput<E>, b: &E::ScalarField, c: &E::G1) {
    let b_final = evaluate_final_poly_for_g1(&proof.challenges.xs, b);
    verify_with_final_exponent(vk, proof, a_comm, &b_final, c)
}

pub fn verify_unstructured<E: Pairing>(vk: &VerifierKey<E>, proof: &Proof<E>, a_comm: &PairingOutput<E>, b: &[E::ScalarField], c: &E::G1) {
    let coeffs = final_folding_exponents(&proof.challenges.xs);
    let b_final = coeffs.into_iter().zip(b).map(|(c, b)| c * b).sum();
    verify_with_final_exponent(vk, proof, a_comm, &b_final, c)
}

fn verify_with_final_exponent<E: Pairing>(vk: &VerifierKey<E>, proof: &Proof<E>, a_comm: &PairingOutput<E>, b_final: &E::ScalarField, c: &E::G1) {
    let challenges = &proof.challenges;

    let xs = challenges.xs.to_vec();
    let mut xs_inv = xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let h1 = (vk.h1 * challenges.c1).into_affine();

    let z = challenges.z;
    let fv_at_z = evaluate_final_poly_for_g2(&xs, &z);
    assert!(kzg::verify_g2(&vk.kzg_vk_g2, proof.v_final, z, fv_at_z, proof.kzg_proof_g2));

    let exps = [xs_inv.as_slice(), xs.as_slice()].concat();
    let bases = [proof.l_comms.as_slice(), proof.r_comms.as_slice()].concat();
    assert_eq!(exps.len(), bases.len());
    let comm = PairingOutput::msm(&bases, &exps).unwrap();
    // TODO: optimize pairings
    let extra = E::pairing(c, h1);
    let comm = comm + a_comm + extra;
    let c_final = proof.a_final * b_final;
    // TODO: batch conversion
    assert_eq!(comm, E::multi_pairing([proof.a_final, c_final.into_affine()],
                                      [proof.v_final, h1]));
}


#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_std::{test_rng, UniformRand};

    use super::*;

    fn _test_mipp_k_for_powers<E: Pairing>() {
        let rng = &mut test_rng();

        let log_n = 8;
        let n = 2usize.pow(log_n);

        let (pk, vk) = setup::<E>(log_n);

        // Want to prove  <A, pow(b)> = A1 + bA2 + ... + b^(n-1)An = C
        let a: Vec<E::G1Affine> = (0..n).map(|_| E::G1Affine::rand(rng)).collect();
        let b = E::ScalarField::rand(rng);
        let c: E::G1 = VariableBaseMSM::msm(&a, &powers(b, n)).unwrap();

        let v: Vec<E::G2Affine> = pk.ck_g2.iter().cloned().step_by(2).collect();

        // A_comm = <A, V> = e(A1, V1) * ... * e(An, Vn)
        let a_comm: PairingOutput<E> = E::multi_pairing(&a, v);

        let t_prove = start_timer!(|| format!("MIPP-k-for-powers, log(n) = {}", log_n));
        let proof = prove_for_powers(&pk, &a, b);
        end_timer!(t_prove);

        verify_for_powers(&vk, &proof, &a_comm, &b, &c);
    }

    fn _test_mipp_k_unstructured<E: Pairing>() {
        let rng = &mut test_rng();

        let log_n = 8;
        let n = 2usize.pow(log_n);

        let (pk, vk) = setup::<E>(log_n);

        // Want to prove  <A, b> = b1A1 + b2bA2 + ... + bnAn = C
        let a: Vec<E::G1Affine> = (0..n).map(|_| E::G1Affine::rand(rng)).collect();
        let b: Vec<E::ScalarField> = (0..n).map(|_| E::ScalarField::rand(rng)).collect();
        let c: E::G1 = VariableBaseMSM::msm(&a, &b).unwrap();

        let v: Vec<E::G2Affine> = pk.ck_g2.iter().cloned().step_by(2).collect();

        // A_comm = <A, V> = e(A1, V1) * ... * e(An, Vn)
        let a_comm: PairingOutput<E> = E::multi_pairing(&a, v);

        let t_prove = start_timer!(|| format!("MIPP-k-unstructured, log(n) = {}", log_n));
        let proof = prove_unstructured(&pk, &a, &b);
        end_timer!(t_prove);

        verify_unstructured(&vk, &proof, &a_comm, &b, &c);
    }

    #[test]
    fn test_mipp_k_for_powers() {
        _test_mipp_k_for_powers::<Bls12_381>();
    }

    #[test]
    fn test_mipp_k_unstructured() {
        _test_mipp_k_for_powers::<Bls12_381>();
    }
}