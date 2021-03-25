use ark_poly::{Radix2EvaluationDomain, EvaluationDomain};
use ark_bw6_761::{BW6_761, Fr};
use ark_ec::{ProjectiveCurve, AffineCurve};
use ark_ff::{One, PrimeField, Field, Zero};
use bench_utils::{end_timer, start_timer};
use merlin::Transcript;

use crate::{endo, Proof, utils, KZG_BW6, point_in_g1_complement, Bitmask};
use crate::transcript::ApkTranscript;
use crate::signer_set::SignerSetCommitment;
use crate::kzg::{VerifierKey, PreparedVerifierKey};
use crate::bls::PublicKey;
use crate::fsrng::fiat_shamir_rng;
use ark_ec::short_weierstrass_jacobian::GroupProjective;
use crate::constraints::{Constraints, SuccinctlyAccountableRegisters, SuccinctAccountableRegisterEvaluations, RegisterEvaluations};


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
        proof: &Proof<SuccinctAccountableRegisterEvaluations>,
    ) -> bool
    {
        assert_eq!(bitmask.size(), self.pks_comm.signer_set_size);

        let mut transcript = self.preprocessed_transcript.clone();
        transcript.append_public_input(&apk, bitmask);
        transcript.append_proof_point(b"b_comm", &proof.register_commitments.basic_commitments.b_comm);
        transcript.append_proof_point(b"acc_x_comm", &proof.register_commitments.basic_commitments.acc_comm.0);
        transcript.append_proof_point(b"acc_y_comm", &proof.register_commitments.basic_commitments.acc_comm.1);
        let r = transcript.get_128_bit_challenge(b"r"); // bitmask batching challenge
        transcript.append_proof_point(b"c_comm", &proof.register_commitments.c_comm);
        transcript.append_proof_point(b"acc_comm", &proof.register_commitments.acc_comm);
        let phi = transcript.get_128_bit_challenge(b"phi"); // constraint polynomials batching challenge
        transcript.append_proof_point(b"q_comm", &proof.q_comm);
        let zeta = transcript.get_128_bit_challenge(b"zeta"); // evaluation point challenge

        let b = proof.register_evaluations.get_bitmask();
        let t_linear_accountability = start_timer!(|| "linear accountability check");
        let b_at_zeta = utils::barycentric_eval_binary_at(zeta, &bitmask, self.domain);
        assert_eq!(b_at_zeta, b);
        // accountability
        end_timer!(t_linear_accountability);

        transcript.append_evals(&proof.register_evaluations);
        transcript.append_proof_scalar(b"q_zeta", &proof.q_zeta);
        transcript.append_proof_scalar(b"r_zeta_omega", &proof.r_zeta_omega);
        let nu: Fr = transcript.get_128_bit_challenge(b"nu"); // KZG opening batching challenge

        let evals_at_zeta = utils::lagrange_evaluations(zeta, self.domain);

        // TODO: 128-bit mul
        let r_comm = proof.register_evaluations.restore_commitment_to_linearization_polynomial(
            phi,
            evals_at_zeta.zeta_minus_omega_inv,
            &proof.register_commitments
        ).into_affine();

        let t_multiexp = start_timer!(|| "multiexp");
        let mut  commitments = vec![
            self.pks_comm.pks_x_comm,
            self.pks_comm.pks_y_comm,
        ];
        commitments.extend(proof.register_commitments.as_vec());
        commitments.push(proof.q_comm);
        let w_comm = KZG_BW6::aggregate_commitments(nu, &commitments);
        end_timer!(t_multiexp);

        let t_opening_points = start_timer!(|| "opening points evaluation");
        let mut register_evals = proof.register_evaluations.as_vec();
        register_evals.push(proof.q_zeta);
        let w_at_zeta = KZG_BW6::aggregate_values(nu, &register_evals);
        end_timer!(t_opening_points);

        let t_kzg_batch_opening = start_timer!(|| "batched KZG openning");
        transcript.append_proof_point(b"w_at_zeta_proof", &proof.w_at_zeta_proof);
        transcript.append_proof_point(b"r_at_zeta_omega_proof", &proof.r_at_zeta_omega_proof);
        let fsrng = &mut fiat_shamir_rng(&mut transcript);
        let (total_c, total_w) = KZG_BW6::aggregate_openings(&self.kzg_pvk,
                                                             &[w_comm, r_comm],
                                                             &[zeta, evals_at_zeta.zeta_omega],
                                                             &[w_at_zeta, proof.r_zeta_omega],
                                                             &[proof.w_at_zeta_proof, proof.r_at_zeta_omega_proof],
                                                             fsrng,
        );
        assert!(KZG_BW6::batch_check_aggregated(&self.kzg_pvk, total_c, total_w));
        end_timer!(t_kzg_batch_opening);

        let t_lazy_subgroup_checks = start_timer!(|| "2 point lazy subgroup check");
        endo::subgroup_check(&total_c);
        endo::subgroup_check(&total_w);
        end_timer!(t_lazy_subgroup_checks);

        let bits_in_bitmask_chunk = 256;
        let bits_in_big_int_limb = 64;
        assert_eq!(bits_in_bitmask_chunk % bits_in_big_int_limb, 0);
        let limbs_in_chunk = bits_in_bitmask_chunk / bits_in_big_int_limb;
        assert_eq!(self.domain.size % bits_in_bitmask_chunk, 0);
        let chunks_in_bitmask = self.domain.size / bits_in_bitmask_chunk; // TODO: bitmask should be right-padded with 0s to domain_size

        let bits_in_bitmask_chunk_inv = Fr::from(256u16).inverse().unwrap();

        let powers_of_r = utils::powers(r, (chunks_in_bitmask - 1) as usize);
        let r_pow_m = r * powers_of_r.last().unwrap();
        let bitmask_chunks = bitmask.to_chunks_as_field_elements::<Fr>(limbs_in_chunk as usize);
        assert_eq!(powers_of_r.len(), bitmask_chunks.len());
        let aggregated_bitmask = bitmask_chunks.into_iter()
            .zip(powers_of_r)
            .map(|(bj, rj)| bj * rj)
            .sum::<Fr>();


        let t_a_zeta_omega1 = start_timer!(|| "A(zw) as fraction");
        let zeta_omega_pow_m = evals_at_zeta.zeta_omega.pow([chunks_in_bitmask]); // m = chunks_in_bitmask
        let zeta_omega_pow_n = zeta_omega_pow_m.pow([bits_in_bitmask_chunk]); // n = domain_size
        let a_zeta_omega1 = bits_in_bitmask_chunk_inv * (zeta_omega_pow_n - Fr::one()) / (zeta_omega_pow_m - Fr::one());
        end_timer!(t_a_zeta_omega1);

        let t_a_zeta_omega2 = start_timer!(|| "A(zw) as polynomial");
        let zeta_omega_pow_m = evals_at_zeta.zeta_omega.pow([chunks_in_bitmask]); // m = chunks_in_bitmask
        let a_zeta_omega2 = bits_in_bitmask_chunk_inv * utils::powers(zeta_omega_pow_m, (bits_in_bitmask_chunk - 1) as usize).iter().sum::<Fr>();
        end_timer!(t_a_zeta_omega2);

        assert_eq!(a_zeta_omega1, a_zeta_omega2);
        let two = Fr::from(2u8);
        let a = two + (r / two.pow([255u64]) - two) * a_zeta_omega1;

        let apk = apk.0.into_affine();

        let constraint_polynomial_evals = proof.register_evaluations.evaluate_constraint_polynomials(apk, &evals_at_zeta, a, r_pow_m, aggregated_bitmask);
        let w = utils::horner_field(&constraint_polynomial_evals, phi);
        proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial
    }
}

