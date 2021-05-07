use ark_poly::{Radix2EvaluationDomain, EvaluationDomain};
use ark_bw6_761::{BW6_761, Fr};
use ark_ec::ProjectiveCurve;
use ark_std::{end_timer, start_timer};
use merlin::Transcript;

use crate::{endo, Proof, utils, KZG_BW6, point_in_g1_complement, Bitmask, RegisterCommitments};
use crate::transcript::ApkTranscript;
use crate::signer_set::SignerSetCommitment;
use crate::kzg::{VerifierKey, PreparedVerifierKey};
use crate::bls::PublicKey;
use crate::fsrng::fiat_shamir_rng;
use crate::piop::bitmask_packing::{SuccinctAccountableRegisterEvaluations, BitmaskPackingCommitments};
use crate::piop::{RegisterPolynomials, RegisterEvaluations};
use crate::piop::affine_addition::{AffineAdditionEvaluations, PartialSumsCommitments};


pub struct Verifier {
    domain: Radix2EvaluationDomain<Fr>,
    kzg_pvk: PreparedVerifierKey<BW6_761>,
    h: ark_bls12_377::G1Affine,
    pks_comm: SignerSetCommitment,
    preprocessed_transcript: Transcript,
}

impl Verifier {

    pub fn verify_simple(
        &self,
        apk: &PublicKey,
        bitmask: &Bitmask,
        proof: &mut Proof<AffineAdditionEvaluations, PartialSumsCommitments, ()>
    ) -> bool {
        self.verify::<
            (),
            PartialSumsCommitments,
            AffineAdditionEvaluations
        >(apk, bitmask, proof)
    }

    pub fn verify_packed(
        &self,
        apk: &PublicKey,
        bitmask: &Bitmask,
        proof: &mut Proof<SuccinctAccountableRegisterEvaluations, PartialSumsCommitments, BitmaskPackingCommitments>
    ) -> bool {
        self.verify::<
            BitmaskPackingCommitments,
            PartialSumsCommitments,
            SuccinctAccountableRegisterEvaluations
        >(apk, bitmask, proof)
    }

    fn verify<AC, C, E>(
        &self,
        apk: &PublicKey,
        bitmask: &Bitmask,
        proof: &mut Proof<E, C, AC>,
    ) -> bool
    where
        AC: RegisterCommitments,
        C: RegisterCommitments,
        E: RegisterEvaluations<C = C> + RegisterEvaluations<AC = AC>,
    {
        assert_eq!(bitmask.size(), self.pks_comm.signer_set_size);

        let mut transcript = self.preprocessed_transcript.clone();

        transcript.append_public_input(&apk, bitmask);
        transcript.append_basic_commitments(&proof.register_commitments);
        let r = transcript.get_128_bit_challenge(b"r"); // bitmask batching challenge
        transcript.append_accountability_commitments(&proof.additional_commitments);
        let phi = transcript.get_128_bit_challenge(b"phi"); // constraint polynomials batching challenge
        transcript.append_proof_point(b"q_comm", &proof.q_comm);
        let zeta = transcript.get_128_bit_challenge(b"zeta"); // evaluation point challenge
        transcript.append_evals(&proof.register_evaluations);
        transcript.append_proof_scalar(b"q_zeta", &proof.q_zeta);
        transcript.append_proof_scalar(b"r_zeta_omega", &proof.r_zeta_omega);
        let nu: Fr = transcript.get_128_bit_challenge(b"nu"); // KZG opening batching challenge

        let evals_at_zeta = utils::lagrange_evaluations(zeta, self.domain);

        proof.register_evaluations.set_bitmask_at_zeta(|| {
            let t_linear_accountability = start_timer!(|| "linear accountability check");
            let b_at_zeta = utils::barycentric_eval_binary_at(zeta, &bitmask, self.domain);
            end_timer!(t_linear_accountability);
            b_at_zeta
        });

        let t_kzg = start_timer!(|| "KZG check");
        // Reconstruct the commitment to the linearization polynomial using the commitments to the registers from the proof.
        let t_r_comm = start_timer!(|| "linearization polynomial commitment");
        // TODO: 128-bit mul
        let r_comm = proof.register_evaluations.restore_commitment_to_linearization_polynomial(
            phi,
            evals_at_zeta.zeta_minus_omega_inv,
            &proof.register_commitments,
            &proof.additional_commitments,
        ).into_affine();
        end_timer!(t_r_comm);


        // Aggregate the commitments to be opened in \zeta, using the challenge \nu.
        let t_multiexp = start_timer!(|| "aggregated commitment");
        let mut  commitments = vec![
            self.pks_comm.pks_x_comm,
            self.pks_comm.pks_y_comm,
        ];
        commitments.extend(proof.register_commitments.as_vec());
        commitments.extend(proof.additional_commitments.as_vec());
        commitments.push(proof.q_comm);
        let w_comm = KZG_BW6::aggregate_commitments(nu, &commitments);
        end_timer!(t_multiexp);


        let t_opening_points = start_timer!(|| "aggregated evaluation");
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


        let t_lazy_subgroup_checks = start_timer!(|| "lazy subgroup check");
        endo::subgroup_check(&total_c);
        endo::subgroup_check(&total_w);
        end_timer!(t_lazy_subgroup_checks);
        end_timer!(t_kzg);

        let apk = apk.0.into_affine();
        let constraint_polynomial_evals = proof.register_evaluations.evaluate_constraint_polynomials(apk, &evals_at_zeta, r, bitmask, self.domain.size);
        let w = utils::horner_field(&constraint_polynomial_evals, phi);
        proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial
    }

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
}

