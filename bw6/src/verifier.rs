use ark_poly::{Radix2EvaluationDomain, EvaluationDomain};
use ark_bw6_761::{BW6_761, Fr};
use ark_ec::ProjectiveCurve;
use ark_ff::{One, PrimeField};
use bench_utils::{end_timer, start_timer};
use merlin::Transcript;

use crate::{endo, Proof, utils, KZG_BW6, point_in_g1_complement, Bitmask};
use crate::transcript::ApkTranscript;
use crate::signer_set::SignerSetCommitment;
use crate::kzg::{VerifierKey, PreparedVerifierKey};
use crate::bls::PublicKey;
use crate::fsrng::fiat_shamir_rng;


pub struct Verifier {
    domain: Radix2EvaluationDomain<Fr>,
    kzg_pvk: PreparedVerifierKey<BW6_761>,
    h: ark_bls12_377::G1Affine,
    pks_comm: SignerSetCommitment,
    preprocessed_transcript: Transcript,
}

impl Verifier {
    pub fn new(
        domain_size: usize,
        kzg_vk: VerifierKey<BW6_761>,
        pks_comm: SignerSetCommitment,
        mut empty_transcript: Transcript,
    ) -> Self {
        // empty_transcript.set_protocol_params(); //TODO
        empty_transcript.set_signer_set(&pks_comm);
        let domain = Radix2EvaluationDomain::<Fr>::new(domain_size).unwrap();
        let kzg_pvk = kzg_vk.prepare();
        Self { domain, kzg_pvk, h: point_in_g1_complement(), pks_comm, preprocessed_transcript: empty_transcript }
    }

    pub fn verify(
        &self,
        apk: &PublicKey,
        bitmask: &Bitmask,
        proof: &Proof,
    ) -> bool
    {
        assert_eq!(bitmask.size(), self.pks_comm.signer_set_size);

        let mut transcript = self.preprocessed_transcript.clone();
        transcript.append_public_input(&apk, bitmask);
        transcript.append_proof_point(b"b_comm", &proof.b_comm);
        transcript.append_proof_point(b"acc_x_comm", &proof.acc_x_comm);
        transcript.append_proof_point(b"acc_y_comm", &proof.acc_y_comm);
        let phi = transcript.get_128_bit_challenge(b"phi"); // constraint polynomials batching challenge
        transcript.append_proof_point(b"q_comm", &proof.q_comm);
        let zeta = transcript.get_128_bit_challenge(b"zeta");

        let t_accountability = start_timer!(|| "accountability check");
        let b_at_zeta = utils::barycentric_eval_binary_at(zeta, &bitmask, self.domain);
        assert_eq!(b_at_zeta, proof.b_zeta);
        // accountability
        end_timer!(t_accountability);

        transcript.append_proof_scalar(b"b_zeta", &proof.b_zeta);
        transcript.append_proof_scalar(b"pks_x_zeta", &proof.pks_x_zeta);
        transcript.append_proof_scalar(b"pks_y_zeta", &proof.pks_y_zeta);
        transcript.append_proof_scalar(b"acc_x_zeta", &proof.acc_x_zeta);
        transcript.append_proof_scalar(b"acc_y_zeta", &proof.acc_y_zeta);
        transcript.append_proof_scalar(b"q_zeta", &proof.q_zeta);
        let nu: Fr = transcript.get_128_bit_challenge(b"nu");

        let t_multiexp = start_timer!(|| "multiexp");
        let nu_repr = nu.into_repr();
        let w2_comm = utils::horner(&[proof.acc_x_comm, proof.acc_y_comm], nu_repr).into_affine();
        let w1_comm = utils::horner(&[self.pks_comm.pks_x_comm, self.pks_comm.pks_y_comm, proof.b_comm, proof.q_comm, w2_comm], nu_repr).into_affine();
        end_timer!(t_multiexp);

        let t_opening_points = start_timer!(|| "opening points evaluation");
        let w1_zeta = utils::horner_field(&[proof.pks_x_zeta, proof.pks_y_zeta, proof.b_zeta, proof.q_zeta, proof.acc_x_zeta, proof.acc_y_zeta], nu);
        let zeta_omega = zeta * self.domain.group_gen;
        let w2_zeta_omega = utils::horner_field(&[proof.acc_x_zeta_omega, proof.acc_y_zeta_omega], nu);
        end_timer!(t_opening_points);

        let t_kzg_batch_opening = start_timer!(|| "batched KZG openning");
        let fsrng = &mut fiat_shamir_rng(&mut transcript);
        let (total_c, total_w) = KZG_BW6::aggregate_openings(&self.kzg_pvk,
                                                             &[w1_comm, w2_comm],
                                                             &[zeta, zeta_omega],
                                                             &[w1_zeta, w2_zeta_omega],
                                                             &[proof.w1_proof, proof.w2_proof],
                                                             fsrng,
        );
        assert!(KZG_BW6::batch_check_aggregated(&self.kzg_pvk, total_c, total_w));
        end_timer!(t_kzg_batch_opening);

        let t_lazy_subgroup_checks = start_timer!(|| "2 point lazy subgroup check");
        endo::subgroup_check(&total_c);
        endo::subgroup_check(&total_w);
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
                ) + (Fr::one() - b) * (y3 - y1);

            let a2 =
                b * (
                    (x1 - x2) * (y3 + y1)
                        - (y2 - y1) * (x3 - x1)
                ) + (Fr::one() - b) * (x3 - x1);

            let a3 = b * (Fr::one() - b);

            let evals = utils::lagrange_evaluations(zeta, self.domain);
            let apk = apk.0.into_affine();
            let apk_plus_h = self.h + apk;
            let a4 = (x1 - self.h.x) * evals.l_0 + (x1 - apk_plus_h.x) * evals.l_minus_1;
            let a5 = (y1 - self.h.y) * evals.l_0 + (y1 - apk_plus_h.y) * evals.l_minus_1;

            let s = zeta - self.domain.group_gen_inv;
            a1 * s + phi * (a2 * s + phi * (a3 + phi * (a4 + phi * a5))) == proof.q_zeta * evals.vanishing_polynomial
        };
    }
}

