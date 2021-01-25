use ark_bw6_761::{BW6_761, Fr as F, G1Projective};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{One, PrimeField};
use ark_std::{UniformRand, test_rng};

use bitvec::vec::BitVec;
use bench_utils::{end_timer, start_timer};

use crate::{endo, Proof, PublicKey, utils, PreparedVerifierKey};
use merlin::Transcript;
use crate::transcript::ApkTranscript;
use crate::signer_set::SignerSetCommitment;

pub struct Verifier {
    vk: PreparedVerifierKey,
    pks_comm: SignerSetCommitment,
    preprocessed_transcript: Transcript,
}

impl Verifier {
    pub fn new(
        vk: PreparedVerifierKey,
        pks_comm: SignerSetCommitment,
        mut empty_transcript: Transcript,
    ) -> Self {
        // empty_transcript.set_protocol_params(); //TODO
        empty_transcript.set_signer_set(&pks_comm);
        Self { vk, pks_comm, preprocessed_transcript: empty_transcript }
    }

    pub fn verify(
        &self,
        apk: &PublicKey,
        bitmask: &BitVec,
        proof: &Proof,
    ) -> bool
    {
        let mut transcript = self.preprocessed_transcript.clone();
        let rng = &mut test_rng(); //TODO: remove
        transcript.append_public_input(&apk, bitmask);
        let f = transcript.get_128_bit_challenge(b"phi");

        transcript.append_proof_point(b"b_comm", &proof.b_comm);
        transcript.append_proof_point(b"acc_x_comm", &proof.acc_x_comm);
        transcript.append_proof_point(b"acc_y_comm", &proof.acc_y_comm);
        transcript.append_proof_point(b"q_comm", &proof.q_comm);
        let zeta = transcript.get_128_bit_challenge(b"zeta");

        let t_accountability = start_timer!(|| "accountability check");
        let b_at_zeta = utils::barycentric_eval_binary_at(zeta, &bitmask, self.vk.domain);
        assert_eq!(b_at_zeta, proof.b_zeta);
        // accountability
        end_timer!(t_accountability);

        transcript.append_proof_scalar(b"b_zeta", &proof.b_zeta);
        transcript.append_proof_scalar(b"pks_x_zeta", &proof.pks_x_zeta);
        transcript.append_proof_scalar(b"pks_y_zeta", &proof.pks_y_zeta);
        transcript.append_proof_scalar(b"acc_x_zeta", &proof.acc_x_zeta);
        transcript.append_proof_scalar(b"acc_y_zeta", &proof.acc_y_zeta);
        transcript.append_proof_scalar(b"q_zeta", &proof.q_zeta);
        let nu: F = transcript.get_128_bit_challenge(b"nu");

        let t_multiexp = start_timer!(|| "multiexp");
        let nu_repr = nu.into_repr();
        let w2_comm = utils::horner(&[proof.acc_x_comm, proof.acc_y_comm], nu_repr).into_affine();
        let w1_comm = utils::horner(&[self.pks_comm.pks_x_comm, self.pks_comm.pks_y_comm, proof.b_comm, proof.q_comm, w2_comm], nu_repr);
        end_timer!(t_multiexp);

        let t_opening_points = start_timer!(|| "opening points evaluation");
        let w1_zeta = utils::horner_field(&[proof.pks_x_zeta, proof.pks_y_zeta, proof.b_zeta, proof.q_zeta, proof.acc_x_zeta, proof.acc_y_zeta], nu);
        let zeta_omega = zeta * self.vk.domain.group_gen;
        let w2_zeta_omega = utils::horner_field(&[proof.acc_x_zeta_omega, proof.acc_y_zeta_omega], nu);
        end_timer!(t_opening_points);

        let r: F = u128::rand(rng).into(); //TODO: deterministic

        let t_kzg_batch_opening = start_timer!(|| "batched KZG openning");
        let c = w1_comm + w2_comm.mul(r); //128-bit mul //TODO: w2_comm is affine
        let v = self.vk.g.mul(w1_zeta + r * w2_zeta_omega); //377-bit FIXED BASE mul
        let z = proof.w1_proof.mul(zeta) + proof.w2_proof.mul(r * zeta_omega); // 128-bit mul + 377 bit mul
        let lhs = c - v + z;
        let mut rhs = proof.w2_proof.mul(r);  //128-bit mul
        rhs.add_assign_mixed(&proof.w1_proof);
        let to_affine = G1Projective::batch_normalization_into_affine(&[lhs, -rhs]); // Basically, not required, BW6 Miller's loop is in projective afair
        let (lhs_affine, rhs_affine) = (to_affine[0], to_affine[1]);
        assert!(BW6_761::product_of_pairings(&[
            (lhs_affine.into(), self.vk.prepared_h.clone()),
            (rhs_affine.into(), self.vk.prepared_beta_h.clone()),
        ]).is_one());
        end_timer!(t_kzg_batch_opening);

        let t_lazy_subgroup_checks = start_timer!(|| "2 point lazy subgroup check");
        endo::subgroup_check(&lhs);
        endo::subgroup_check(&rhs);
        end_timer!(t_lazy_subgroup_checks);

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

            let evals = &self.vk.lagrange_evaluations(zeta);
            let apk = apk.0.into_affine();
            let apk_plus_h = self.vk.h + apk;
            let a4 = (x1 - self.vk.h.x) * evals.l_0 + (x1 - apk_plus_h.x) * evals.l_minus_1;
            let a5 = (y1 - self.vk.h.y) * evals.l_0 + (y1 - apk_plus_h.y) * evals.l_minus_1;

            let s = zeta - self.vk.domain.group_gen_inv;
            a1 * s + f * (a2 * s + f * (a3 + f * (a4 + f * a5))) == proof.q_zeta * evals.vanishing_polynomial
        };
    }
}