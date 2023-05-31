use ark_ec::CurveGroup;
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::VariableBaseMSM;
use ark_ff::{batch_inversion, Field};
use ark_std::{end_timer, start_timer, test_rng, UniformRand};
use ark_std::rand::Rng;

use crate::ipa::{compute_final_poly_for_g1, compute_final_poly_for_g2, evaluate_final_poly_for_g1,
                 evaluate_final_poly_for_g2, fold_points, fold_scalars};
use crate::kzg;


// Inner product argument for <A, b> = b1A1 + ... + bnAn = C, given the commitments:
// A_comm = [A, V] = e(A1, V1) * ... * e(An, Vn)
// b_comm = <b, W> = b1W1 + ... + bnWn, where
// Vi = s^(2i-2)H in G2, and
// Wi = t^(i-1)G in G1.

// { (A_comm, b_comm, C; A, b) | C = <A, b>, A_comm = [A, V], b_comm = <b, W> }


pub struct ProverKey<E: Pairing> {
    log_n: u32,
    // n points in G1: G, tG, ..., t^(n-1)G
    pub ck_g1: Vec<E::G1Affine>,
    // 2n-1 points in G2: H, sH, ..., s^(2n-2)H
    pub ck_g2: Vec<E::G2Affine>,
    h1: E::G2Affine,
    h2: E::G2Affine,
}

pub struct VerifierKey<E: Pairing> {
    kzg_vk_g1: kzg::VerifierKeyG1<E>,
    kzg_vk_g2: kzg::VerifierKeyG2<E>,
    h1: E::G2Affine,
    h2: E::G2Affine,
}

pub struct Proof<E: Pairing> {
    // Left cross-commitments
    l_comms: Vec<PairingOutput<E>>,
    // Right cross-commitments
    r_comms: Vec<PairingOutput<E>>,
    // Final folded point
    a_final: E::G1Affine,
    // Final folded scalar
    b_final: E::ScalarField,
    // Final commitment key in G2
    v_final: E::G2Affine,
    // Final commitment key in G1
    w_final: E::G1Affine,
    // Final commitment key argument proof for the key in G1
    kzg_proof_g1: E::G1Affine,
    // Final commitment key argument proof for the key in G2
    kzg_proof_g2: E::G2Affine,

    challenges: Challenges<E>,
}

struct Challenges<E: Pairing> {
    xs: Vec<E::ScalarField>,
    c1: E::ScalarField,
    c2: E::ScalarField,
    z: E::ScalarField,
}

impl<E: Pairing> Challenges<E> {
    fn new(log_n: usize) -> Self {
        let rng = &mut test_rng();
        Challenges {
            xs: (0..log_n).map(|_| E::ScalarField::from(rng.gen::<u128>())).collect(),
            c1: E::ScalarField::rand(rng),
            c2: E::ScalarField::rand(rng),
            z: E::ScalarField::rand(rng),
        }
    }
}

pub fn setup<E: Pairing>(log_n: u32) -> (ProverKey<E>, VerifierKey<E>) {
    let rng = &mut test_rng();

    let n = 2usize.pow(log_n);
    let g = E::G1::rand(rng);
    let h = E::G2::rand(rng);

    let (ck_g1, kzg_vk_g1) = kzg::setup_g1::<E>(n, g, h);
    let (ck_g2, kzg_vk_g2) = kzg::setup_g2::<E>(2 * n - 1, g, h);

    let h1 = E::G2Affine::rand(rng);
    let h2 = E::G2Affine::rand(rng);

    let pk = ProverKey {
        log_n,
        ck_g1,
        ck_g2,
        h1,
        h2,
    };

    let vk = VerifierKey {
        kzg_vk_g1,
        kzg_vk_g2,
        h1,
        h2,
    };

    (pk, vk)
}

pub fn prove<E: Pairing>(pk: &ProverKey<E>, a: &[E::G1Affine], b: &[E::ScalarField]) -> Proof<E> {
    let log_n = pk.log_n as usize;
    let n = 2usize.pow(log_n as u32);

    // TODO: Fiat-Shamir
    let challenges = Challenges::<E>::new(log_n);

    let h1 = (pk.h1 * challenges.c1).into_affine();
    let h2 = (pk.h2 * challenges.c2).into_affine();

    let mut n1 = n;
    let mut a_folded = a.to_vec();
    let mut b_folded = b.to_vec();
    let mut v_folded: Vec<E::G2Affine> = pk.ck_g2.iter().cloned().step_by(2).collect();
    let mut w_folded: Vec<E::G1Affine> = pk.ck_g1.clone();
    assert_eq!(a_folded.len(), n);
    assert_eq!(b_folded.len(), n);
    assert_eq!(v_folded.len(), n);
    assert_eq!(w_folded.len(), n);

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
        let wl = &w_folded[..n1];
        let wr = &w_folded[n1..];

        let t_msm = start_timer!(|| format!("4 x {}-msm in G1", n1));
        let bl_comm = E::G1::msm(wr, bl).unwrap();
        let br_comm = E::G1::msm(wl, br).unwrap();
        let cl = E::G1::msm(ar, bl).unwrap();
        let cr = E::G1::msm(al, br).unwrap();
        end_timer!(t_msm);

        // TODO: batch conversion to affine
        let bl_comm = bl_comm.into_affine();
        let br_comm = br_comm.into_affine();
        let cl = cl.into_affine();
        let cr = cr.into_affine();

        let l_vals = [ar, &[bl_comm, cl]].concat();
        let r_vals = [al, &[br_comm, cr]].concat();
        let l_leys = [vl, &[h1, h2]].concat();
        let r_keys = [vr, &[h1, h2]].concat();

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
        w_folded = fold_points(wl, wr, &x_inv);
        end_timer!(t_folding);

        end_timer!(t_round);
    }
    end_timer!(t_gipa);

    assert_eq!(a_folded.len(), 1);
    assert_eq!(b_folded.len(), 1);
    assert_eq!(v_folded.len(), 1);
    assert_eq!(w_folded.len(), 1);
    let a_final = a_folded[0].into();
    let b_final = b_folded[0];
    let v_final = v_folded[0];
    let w_final = w_folded[0];

    let mut xs_inv = challenges.xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let f_w = compute_final_poly_for_g1(&xs_inv);
    let f_v = compute_final_poly_for_g2(&challenges.xs);
    let t_kzg = start_timer!(|| format!("{}-msm in G1 + {}-msm in G2", n, 2 * n));
    let kzg_proof_g1 = kzg::open_g1::<E>(&pk.ck_g1, &f_w, challenges.z);
    let kzg_proof_g2 = kzg::open_g2::<E>(&pk.ck_g2, &f_v, challenges.z);
    end_timer!(t_kzg);

    Proof {
        l_comms,
        r_comms,
        a_final,
        b_final,
        v_final,
        w_final,
        kzg_proof_g1,
        kzg_proof_g2,
        challenges,
    }
}

pub fn verify<E: Pairing>(vk: &VerifierKey<E>, proof: &Proof<E>, a_comm: &PairingOutput<E>, b_comm: &E::G1, c: &E::G1) {
    let challenges = &proof.challenges;

    let xs = challenges.xs.to_vec();
    let mut xs_inv = xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let h1 = (vk.h1 * challenges.c1).into_affine();
    let h2 = (vk.h2 * challenges.c2).into_affine();

    let z = challenges.z;
    let fv_at_z = evaluate_final_poly_for_g2(&xs, &z);
    let fw_at_z = evaluate_final_poly_for_g1(&xs_inv, &z);
    assert!(kzg::verify_g2(&vk.kzg_vk_g2, proof.v_final, z, fv_at_z, proof.kzg_proof_g2));
    assert!(kzg::verify_g1(&vk.kzg_vk_g1, proof.w_final, z, fw_at_z, proof.kzg_proof_g1));

    let exps = [xs_inv, xs].concat();
    let bases = [proof.l_comms.as_slice(), proof.r_comms.as_slice()].concat();
    assert_eq!(exps.len(), bases.len());
    let comm = PairingOutput::msm(&bases, &exps).unwrap();
    // TODO: optimize pairings
    let extra = E::multi_pairing([b_comm, c], [h1, h2]);
    let comm = comm + a_comm + extra;
    let b_comm_final = proof.w_final * proof.b_final;
    let c_final = proof.a_final * proof.b_final;
    // TODO: batch conversion
    assert_eq!(comm, E::multi_pairing([proof.a_final, b_comm_final.into_affine(), c_final.into_affine()],
                                      [proof.v_final, h1, h2]));
}


#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_std::{test_rng, UniformRand};

    use super::*;

    fn _test_mipp<E: Pairing>() {
        let rng = &mut test_rng();

        let log_n = 8;
        let n = 2usize.pow(log_n);

        let (pk, vk) = setup::<E>(log_n);

        // Want to prove <A, b> = b1A1 + ... + bnAn = C
        let a: Vec<E::G1Affine> = (0..n).map(|_| E::G1Affine::rand(rng)).collect();
        let b: Vec<E::ScalarField> = (0..n).map(|_| E::ScalarField::rand(rng)).collect();
        let c: E::G1 = VariableBaseMSM::msm(&a, &b).unwrap();

        let v: Vec<E::G2Affine> = pk.ck_g2.iter().cloned().step_by(2).collect();
        let w: Vec<E::G1Affine> = pk.ck_g1.clone();

        // A_comm = <A, V> = e(A1, V1) * ... * e(An, Vn)
        let a_comm: PairingOutput<E> = E::multi_pairing(&a, v);
        // b_comm = <b, W> = b1W1 * ... + bnWn
        let b_comm: E::G1 = VariableBaseMSM::msm(&w, &b).unwrap();

        let t_prove = start_timer!(|| format!("MIPP-u, log(n) = {}", log_n));
        let proof = prove(&pk, &a, &b);
        end_timer!(t_prove);

        verify(&vk, &proof, &a_comm, &b_comm, &c);
    }

    #[test]
    fn test_mipp() {
        _test_mipp::<Bls12_381>();
    }
}