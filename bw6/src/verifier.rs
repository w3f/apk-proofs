use ark_poly::Radix2EvaluationDomain;
use ark_bw6_761::{BW6_761, Fr};
use ark_ec::{ProjectiveCurve, AffineCurve};
use ark_std::{end_timer, start_timer};
use merlin::{Transcript, TranscriptRng};

use crate::{endo, Proof, utils, RegisterCommitments, PublicInput, AccountablePublicInput, CountingPublicInput, SimpleProof, PackedProof, CountingProof, KeysetCommitment, NewKzgBw6};
use crate::transcript::ApkTranscript;
use crate::fsrng::fiat_shamir_rng;
use crate::piop::bitmask_packing::{SuccinctAccountableRegisterEvaluations, BitmaskPackingCommitments};
use crate::piop::{VerifierProtocol, RegisterEvaluations};
use crate::piop::affine_addition::{AffineAdditionEvaluations, PartialSumsCommitments, PartialSumsAndBitmaskCommitments};
use crate::piop::basic::AffineAdditionEvaluationsWithoutBitmask;
use crate::utils::LagrangeEvaluations;
use crate::piop::counting::{CountingEvaluations, CountingCommitments};
use fflonk::aggregation::single::aggregate_claims_multiexp;
use fflonk::pcs::kzg::KzgOpening;
use fflonk::pcs::kzg::params::{KzgVerifierKey, RawKzgVerifierKey};
use ark_ff::{UniformRand, One};


pub struct Verifier {
    domain: Radix2EvaluationDomain<Fr>,
    kzg_pvk: KzgVerifierKey<BW6_761>,
    pks_comm: KeysetCommitment,
    preprocessed_transcript: Transcript,
}

struct Challenges {
    r: Fr,
    phi: Fr,
    zeta: Fr,
    nus: Vec<Fr>,
}


impl Verifier {
    pub fn verify_simple(
        &self,
        public_input: &AccountablePublicInput,
        proof: &SimpleProof,
    ) -> bool {
        assert_eq!(public_input.bitmask.size(), self.pks_comm.keyset_size);
        let (challenges, mut fsrng) = self.restore_challenges(public_input, proof, AffineAdditionEvaluations::POLYS_OPENED_AT_ZETA);
        let evals_at_zeta = utils::lagrange_evaluations(challenges.zeta, self.domain);

        let t_linear_accountability = start_timer!(|| "linear accountability check");
        let b_at_zeta = utils::barycentric_eval_binary_at(challenges.zeta, &public_input.bitmask, self.domain);
        end_timer!(t_linear_accountability);

        let evaluations_with_bitmask = AffineAdditionEvaluations {
            keyset: proof.register_evaluations.keyset,
            bitmask: b_at_zeta,
            partial_sums: proof.register_evaluations.partial_sums,
        };

        self.validate_evaluations::<
            (),
            PartialSumsCommitments,
            AffineAdditionEvaluationsWithoutBitmask,
            AffineAdditionEvaluations,
        >(proof, &evaluations_with_bitmask, &challenges, &mut fsrng, &evals_at_zeta);

        let apk = public_input.apk;
        let constraint_polynomial_evals = evaluations_with_bitmask.evaluate_constraint_polynomials(apk, &evals_at_zeta);
        let w = utils::horner_field(&constraint_polynomial_evals, challenges.phi);
        proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial
    }

    pub fn verify_packed(
        &self,
        public_input: &AccountablePublicInput,
        proof: &PackedProof,
    ) -> bool {
        assert_eq!(public_input.bitmask.size(), self.pks_comm.keyset_size);
        let (challenges, mut fsrng) = self.restore_challenges(public_input, proof, SuccinctAccountableRegisterEvaluations::POLYS_OPENED_AT_ZETA);
        let evals_at_zeta = utils::lagrange_evaluations(challenges.zeta, self.domain);

        self.validate_evaluations::<
            BitmaskPackingCommitments,
            PartialSumsAndBitmaskCommitments,
            SuccinctAccountableRegisterEvaluations,
            SuccinctAccountableRegisterEvaluations,
        >(proof, &proof.register_evaluations, &challenges, &mut fsrng, &evals_at_zeta);

        let apk = public_input.apk;
        let constraint_polynomial_evals = proof.register_evaluations.evaluate_constraint_polynomials(apk, &evals_at_zeta, challenges.r, &public_input.bitmask, self.domain.size);
        let w = utils::horner_field(&constraint_polynomial_evals, challenges.phi);
        proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial
    }

    pub fn verify_counting(
        &self,
        public_input: &CountingPublicInput,
        proof: &CountingProof,
    ) -> bool {
        assert!(public_input.count > 0);
        let (challenges, mut fsrng) = self.restore_challenges(public_input, proof, CountingEvaluations::POLYS_OPENED_AT_ZETA);
        let evals_at_zeta = utils::lagrange_evaluations(challenges.zeta, self.domain);
        let count = Fr::from(public_input.count as u16);

        self.validate_evaluations::<
            (),
            CountingCommitments,
            CountingEvaluations,
            CountingEvaluations,
        >(proof, &proof.register_evaluations, &challenges, &mut fsrng, &evals_at_zeta);

        let apk = public_input.apk;
        let constraint_polynomial_evals = proof.register_evaluations.evaluate_constraint_polynomials(apk, count, &evals_at_zeta);
        let w = utils::horner_field(&constraint_polynomial_evals, challenges.phi);
        proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial
    }


    fn validate_evaluations<AC, C, E, P>(
        &self,
        proof: &Proof<E, C, AC>,
        protocol: &P,
        challenges: &Challenges,
        fsrng: &mut TranscriptRng,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
    ) -> ()
        where
            AC: RegisterCommitments,
            C: RegisterCommitments,
            E: RegisterEvaluations,
            P: VerifierProtocol<C1=C> + VerifierProtocol<C2=AC>,
    {
        let t_kzg = start_timer!(|| "KZG check");
        // Reconstruct the commitment to the linearization polynomial using the commitments to the registers from the proof.
        let t_r_comm = start_timer!(|| "linearization polynomial commitment");
        // TODO: 128-bit mul
        let r_comm = protocol.restore_commitment_to_linearization_polynomial(
            challenges.phi,
            evals_at_zeta.zeta_minus_omega_inv,
            &proof.register_commitments,
            &proof.additional_commitments,
        ).into_affine();
        end_timer!(t_r_comm);


        // Aggregate the commitments to be opened in \zeta, using the challenge \nu.
        let t_aggregate_claims = start_timer!(|| "aggregate evaluation claims in zeta");
        let mut commitments = vec![
            self.pks_comm.pks_comm.0,
            self.pks_comm.pks_comm.1,
        ];
        commitments.extend(proof.register_commitments.as_vec());
        commitments.extend(proof.additional_commitments.as_vec());
        commitments.push(proof.q_comm);
        // ...together with the corresponding values
        let mut register_evals = proof.register_evaluations.as_vec();
        register_evals.push(proof.q_zeta);
        assert_eq!(commitments.len(), challenges.nus.len());
        assert_eq!(register_evals.len(), challenges.nus.len());
        let (w_comm, w_at_zeta) = aggregate_claims_multiexp(commitments, register_evals, &challenges.nus);
        end_timer!(t_aggregate_claims);

        let t_kzg_batch_opening = start_timer!(|| "batched KZG openning");
        let opening_at_zeta = KzgOpening {
            c: w_comm,
            x: challenges.zeta,
            y: w_at_zeta,
            proof: proof.w_at_zeta_proof,
        };
        let opening_at_zeta_omega = KzgOpening {
            c: r_comm,
            x: evals_at_zeta.zeta_omega,
            y: proof.r_zeta_omega,
            proof: proof.r_at_zeta_omega_proof,
        };
        let openings = vec![opening_at_zeta, opening_at_zeta_omega];
        let coeffs = [Fr::one(), u128::rand(fsrng).into()];
        let acc_opening = NewKzgBw6::accumulate(openings, &coeffs, &self.kzg_pvk);
        assert!(NewKzgBw6::verify_accumulated(acc_opening.clone(), &self.kzg_pvk), "KZG verification");
        end_timer!(t_kzg_batch_opening);

        let t_lazy_subgroup_checks = start_timer!(|| "lazy subgroup check");
        assert!(endo::subgroup_check(&acc_opening.acc.into_projective()));
        assert!(endo::subgroup_check(&acc_opening.proof.into_projective()));
        end_timer!(t_lazy_subgroup_checks);

        end_timer!(t_kzg);
    }

    fn restore_challenges<E, C, AC>(&self, public_input: &impl PublicInput, proof: &Proof<E, C, AC>, batch_size: usize) -> (Challenges, TranscriptRng)
        where
            AC: RegisterCommitments,
            C: RegisterCommitments,
            E: RegisterEvaluations,
    {
        let mut transcript = self.preprocessed_transcript.clone();
        transcript.append_public_input(public_input);
        transcript.append_register_commitments(&proof.register_commitments);
        let r = transcript.get_bitmask_aggregation_challenge();
        transcript.append_2nd_round_register_commitments(&proof.additional_commitments);
        let phi = transcript.get_constraints_aggregation_challenge();
        transcript.append_quotient_commitment(&proof.q_comm);
        let zeta = transcript.get_evaluation_point();
        transcript.append_evaluations(&proof.register_evaluations, &proof.q_zeta, &proof.r_zeta_omega);
        let nus = transcript.get_kzg_aggregation_challenges(batch_size);
        (Challenges { r, phi, zeta, nus }, fiat_shamir_rng(&mut transcript))
    }

    pub fn new(
        kzg_vk: RawKzgVerifierKey<BW6_761>,
        pks_comm: KeysetCommitment,
        mut empty_transcript: Transcript,
    ) -> Self {
        empty_transcript.set_protocol_params(&pks_comm.domain, &kzg_vk);
        empty_transcript.set_keyset_commitment(&pks_comm);

        let kzg_pvk = kzg_vk.prepare();
        Self {
            domain: pks_comm.domain,
            kzg_pvk,
            pks_comm,
            preprocessed_transcript: empty_transcript,
        }
    }
}

