use ark_bw6_761::BW6_761;
use ark_poly::{Polynomial, EvaluationDomain};
use merlin::Transcript;

use crate::{KZG_BW6, Proof, Bitmask, PublicInput, SimpleProof, PackedProof, CountingProof, AccountablePublicInput, CountingPublicInput, KeysetCommitment, kzg};
use crate::transcript::ApkTranscript;
use crate::kzg::ProverKey;
use crate::piop::ProverProtocol;
use crate::piop::RegisterPolynomials;
use crate::piop::packed::PackedRegisterBuilder;
use crate::piop::basic::BasicRegisterBuilder;
use crate::piop::counting::CountingScheme;
use crate::keyset::Keyset;
use ark_ec::ProjectiveCurve;


pub struct Prover {
    keyset: Keyset,
    kzg_pk: ProverKey<BW6_761>,
    preprocessed_transcript: Transcript,
}


impl Prover {

    pub fn new(
        mut keyset: Keyset,
        keyset_comm: &KeysetCommitment,
        // prover needs both KZG pk and vk, as it commits to the latter to bind the srs
        kzg_params: kzg::Params<BW6_761>,
        mut empty_transcript: Transcript,
    ) -> Self {
        assert!(kzg_params.fits(keyset.domain.size())); // SRS contains enough elements
        empty_transcript.set_protocol_params(&keyset.domain, &kzg_params.get_vk());
        empty_transcript.set_keyset_commitment(&keyset_comm);

        keyset.amplify();

        Self {
            keyset,
            kzg_pk: kzg_params.get_pk(),
            preprocessed_transcript: empty_transcript,
        }
    }

    pub fn prove_simple(&self, bitmask: Bitmask) -> (SimpleProof, AccountablePublicInput) {
        self.prove::<BasicRegisterBuilder>(bitmask)
    }

    pub fn prove_packed(&self, bitmask: Bitmask) -> (PackedProof, AccountablePublicInput) {
        self.prove::<PackedRegisterBuilder>(bitmask)
    }

    pub fn prove_counting(&self, bitmask: Bitmask) -> (CountingProof, CountingPublicInput) {
        self.prove::<CountingScheme>(bitmask)
    }

    #[allow(non_snake_case)]
    fn prove<P: ProverProtocol>(&self, bitmask: Bitmask) -> (Proof<P::E, <P::P1 as RegisterPolynomials>::C, <P::P2 as RegisterPolynomials>::C>, P::PI)
    {
        assert_eq!(bitmask.size(), self.keyset.size());
        assert!(bitmask.count_ones() > 0); // as EC identity doesn't have and affine representation

        let apk = self.keyset.aggregate(&bitmask.to_bits()).into_affine();

        let mut transcript = self.preprocessed_transcript.clone();
        let public_input = P::PI::new(&apk, &bitmask);
        transcript.append_public_input(&public_input);

        // 1. Compute and commit to the basic registers.
        let mut protocol = P::init(bitmask, self.keyset.clone());
        let partial_sums_polynomials = protocol.get_register_polynomials_to_commit1();
        let partial_sums_commitments = partial_sums_polynomials.commit(
            |p| KZG_BW6::commit(&self.kzg_pk, &p)
        );

        transcript.append_register_commitments(&partial_sums_commitments);

        // 2. Receive bitmask aggregation challenge,
        // compute and commit to succinct accountability registers.
        let r = transcript.get_bitmask_aggregation_challenge();
        // let acc_registers = D::wrap(registers, b, r);
        let acc_register_polynomials = protocol.get_register_polynomials_to_commit2(r);
        let acc_register_commitments = acc_register_polynomials.commit(
            |p| KZG_BW6::commit(&self.kzg_pk, &p)
        );
        transcript.append_2nd_round_register_commitments(&acc_register_commitments);

        // 3. Receive constraint aggregation challenge,
        // compute and commit to the quotient polynomial.
        let phi = transcript.get_constraints_aggregation_challenge();
        let q_poly = protocol.compute_quotient_polynomial(phi, self.keyset.domain);
        let q_comm = KZG_BW6::commit(&self.kzg_pk, &q_poly);
        transcript.append_quotient_commitment(&q_comm);

        // 4. Receive the evaluation point,
        // evaluate register polynomials and the quotient polynomial,
        // compute the linearization polynomial and evaluate it at the shifted evaluation point,
        // commit to all the evaluations.
        let zeta = transcript.get_evaluation_point();
        let register_evaluations = protocol.evaluate_register_polynomials(zeta);
        let q_zeta = q_poly.evaluate(&zeta);
        let zeta_omega = zeta * self.keyset.domain.group_gen;
        let r_poly = protocol.compute_linearization_polynomial(phi, zeta);
        let r_zeta_omega = r_poly.evaluate(&zeta_omega);
        transcript.append_evaluations(&register_evaluations, &q_zeta, &r_zeta_omega);

        // 5. Receive the polynomials aggregation challenge,
        // open the aggregated polynomial at the evaluation point,
        // and the linearization polynomial at the shifted evaluation point,
        // and commit to the opening proofs.
        let nu = transcript.get_kzg_aggregation_challenge();
        let mut register_polynomials = protocol.get_register_polynomials_to_open();
        register_polynomials.push(q_poly);
        let w_poly = KZG_BW6::aggregate_polynomials(nu, &register_polynomials);
        let w_at_zeta_proof = KZG_BW6::open(&self.kzg_pk, &w_poly, zeta);
        let r_at_zeta_omega_proof = KZG_BW6::open(&self.kzg_pk, &r_poly, zeta_omega);

        // Finally, compose the proof.
        let proof = Proof {
            register_commitments: partial_sums_commitments,
            additional_commitments: acc_register_commitments,
            // phi <-
            q_comm,
            // zeta <-
            register_evaluations,
            q_zeta,
            r_zeta_omega,
            // <- nu
            w_at_zeta_proof,
            r_at_zeta_omega_proof,
        };

        (proof, public_input)
    }
}