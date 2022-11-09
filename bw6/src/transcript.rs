use ark_bw6_761::{BW6_761, Fr, G1Affine};
use ark_ff::Field;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use fflonk::pcs::kzg::params::RawKzgVerifierKey;
use merlin::Transcript;

use crate::{KeysetCommitment, PublicInput};
use crate::piop::{RegisterCommitments, RegisterEvaluations};

pub(crate) trait ApkTranscript {

    fn set_protocol_params(&mut self, domain: &Radix2EvaluationDomain<Fr>, kzg_vk: &RawKzgVerifierKey<BW6_761>) {
        self._append_serializable(b"domain", domain);
        self._append_serializable(b"vk", kzg_vk);
    }

    fn set_keyset_commitment(&mut self, keyset_commitment: &KeysetCommitment) {
        self._append_serializable(b"keyset_commitment", keyset_commitment);
    }

    fn append_public_input(&mut self, public_input: &impl PublicInput) {
        self._append_serializable(b"public_input", public_input);
    }

    fn append_register_commitments(&mut self, register_commitments: &impl RegisterCommitments) {
        self._append_serializable(b"register_commitments", register_commitments);
    }

    fn get_bitmask_aggregation_challenge(&mut self) -> Fr {
        self._get_128_bit_challenge(b"bitmask_aggregation")
    }

    fn append_2nd_round_register_commitments(&mut self, register_commitments: &impl RegisterCommitments) {
        self._append_serializable(b"2nd_round_register_commitments", register_commitments);
    }

    fn get_constraints_aggregation_challenge(&mut self) -> Fr {
        self._get_128_bit_challenge(b"constraints_aggregation")
    }

    fn append_quotient_commitment(&mut self, point: &G1Affine) {
        self._append_serializable(b"quotient", point);
    }

    fn get_evaluation_point(&mut self) -> Fr {
        self._get_128_bit_challenge(b"evaluation_point")
    }

    fn append_evaluations(&mut self, evals: &impl RegisterEvaluations, q_at_zeta: &Fr, r_at_zeta_omega: &Fr) {
        self._append_serializable(b"register_evaluations", evals);
        self._append_serializable(b"quotient_evaluation", q_at_zeta);
        self._append_serializable(b"shifted_linearization_evaluation", r_at_zeta_omega);
    }

    fn get_kzg_aggregation_challenges(&mut self, n: usize) -> Vec<Fr> {
        self._get_128_bit_challenges(b"kzg_aggregation", n)
    }

    fn _get_128_bit_challenge(&mut self, label: &'static [u8]) -> Fr;

    fn _get_128_bit_challenges(&mut self, label: &'static [u8], n: usize) -> Vec<Fr>;

    fn _append_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize);
}

impl ApkTranscript for Transcript {

    fn _get_128_bit_challenge(&mut self, label: &'static [u8]) -> Fr {
        let mut buf = [0u8; 16];
        self.challenge_bytes(label, &mut buf);
        Fr::from_random_bytes(&buf).unwrap()
    }

    fn _get_128_bit_challenges(&mut self, label: &'static [u8], n: usize) -> Vec<Fr> {
        (0..n).map(|_| self._get_128_bit_challenge(label)).collect() //TODO: unlikely secure
    }

    fn _append_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize) {
        let mut buf = vec![0; message.compressed_size()];
        message.serialize_compressed(&mut buf).unwrap();
        self.append_message(label, &buf);
    }
}