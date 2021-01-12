use std::time::Instant;

use ark_bw6_761::{BW6_761, Fr as F, G1Projective};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{One, PrimeField, test_rng, UniformRand};

use bitvec::vec::BitVec;
use bench_utils::{end_timer, start_timer};

use crate::{endo, Proof, PublicKey, utils, VerifierKey};

pub fn verify(
    pks_x_comm: &ark_bw6_761::G1Affine,
    pks_y_comm: &ark_bw6_761::G1Affine,
    apk: &PublicKey,
    bitmask: &BitVec,
    proof: &Proof,
    vk: &VerifierKey) -> bool
{
    let rng = &mut test_rng(); //TODO: method parameter

    let nu= proof.nu;

    // let timer = Instant::now();
    let accountability = start_timer!(|| "accountability check");
    let b_at_zeta = utils::barycentric_eval_binary_at(proof.zeta, &bitmask, vk.domain);
    assert_eq!(b_at_zeta, proof.b_zeta); // accountability
    end_timer!(accountability);
    // println!("  {}μs = accountability", timer.elapsed().as_micros());


    let timer = Instant::now();
    let nu_repr = nu.into_repr();
    let w2_comm = utils::horner(&[proof.acc_x_comm, proof.acc_y_comm], nu_repr).into_affine();
    let w1_comm = utils::horner(&[*pks_x_comm, *pks_y_comm, proof.b_comm, proof.q_comm, w2_comm], nu_repr);
    // println!("  {}μs = multiexp", timer.elapsed().as_micros());

    let timer = Instant::now();
    let w1_zeta = utils::horner_field(&[proof.pks_x_zeta, proof.pks_y_zeta, proof.b_zeta, proof.q_zeta, proof.acc_x_zeta, proof.acc_y_zeta], nu);
    let zeta_omega = proof.zeta * vk.domain.group_gen;
    let w2_zeta_omega = utils::horner_field(&[proof.acc_x_zeta_omega, proof.acc_y_zeta_omega], nu);
    // println!("  {}μs = opening points evaluation", timer.elapsed().as_micros());

    let r: F = u128::rand(rng).into();

    let timer = Instant::now();

    let c = w1_comm + w2_comm.mul(r); //128-bit mul //TODO: w2_comm is affine
    let v = vk.g.mul(w1_zeta + r * w2_zeta_omega); //377-bit FIXED BASE mul
    let z = proof.w1_proof.mul(proof.zeta) + proof.w2_proof.mul(r * zeta_omega); // 128-bit mul + 377 bit mul
    let lhs = c - v + z;

    let mut rhs = proof.w2_proof.mul(r);  //128-bit mul
    rhs.add_assign_mixed(&proof.w1_proof);

    let to_affine = G1Projective::batch_normalization_into_affine(&[lhs, -rhs]); // Basically, not required, BW6 Miller's loop is in projective afair
    let (lhs_affine, rhs_affine) = (to_affine[0], to_affine[1]);
    assert!(BW6_761::product_of_pairings(&[
        (lhs_affine.into(), vk.prepared_h.clone()),
        (rhs_affine.into(), vk.prepared_beta_h.clone()),
    ]).is_one());
    // println!("  {}μs = batched KZG openning", timer.elapsed().as_micros());

    let timer = Instant::now();
    endo::subgroup_check(&lhs);
    endo::subgroup_check(&rhs);
    // println!("  {}μs = 2-point subgroup check", timer.elapsed().as_micros());

    return {
        let b = proof.b_zeta;
        let x1 = proof.acc_x_zeta;
        let y1 = proof.acc_y_zeta;
        let x2 = proof.pks_x_zeta;
        let y2 = proof.pks_y_zeta;
        let x3 = proof.acc_x_zeta_omega;
        let y3 = proof.acc_y_zeta_omega;

        let a1 =
            b * (
                (x1 - x2) * (x1 - x2) * (x1 + x2 + x3)
                    - (y2 - y1) * (y2 - y1)
            ) + (F::one() - b) * (y3 - y1);

        let a2 =
            b * (
                (x1 - x2) * (y3 + y1)
                    - (y2 - y1) * (x3 - x1)
            ) + (F::one() - b) * (x3 - x1);

        let a3 = b * (F::one() - b);

        let evals = &vk.lagrange_evaluations(proof.zeta);
        let apk = apk.0.into_affine();
        let apk_plus_h = vk.h + apk;
        let a4 = (x1 - vk.h.x) * evals.l_0 + (x1 - apk_plus_h.x) * evals.l_minus_1;
        let a5 = (y1 - vk.h.y) * evals.l_0 + (y1 - apk_plus_h.y) * evals.l_minus_1;

        let s = proof.zeta - vk.domain.group_gen_inv;
        let f = proof.phi;
        a1 * s + f * (a2 * s + f * (a3 + f * (a4 + f * a5))) == proof.q_zeta * evals.vanishing_polynomial
    }
}