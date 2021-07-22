use merlin::Transcript;
use ark_serialize::CanonicalSerialize;
use ark_ff::{Field, ToBytes};
use ark_ec::ProjectiveCurve;
use crate::signer_set::SignerSetCommitment;
use crate::bls::PublicKey;
use crate::{Bitmask, PublicInput};
use crate::piop::{RegisterCommitments, RegisterEvaluations};
use ark_poly::Radix2EvaluationDomain;
use ark_bw6_761::{Fr, BW6_761};
use crate::kzg::VerifierKey;


pub(crate) trait ApkTranscript {

    fn set_protocol_params(&mut self, domain: &Radix2EvaluationDomain<Fr>, kzg_vk: &VerifierKey<BW6_761>) {
        self._append_serializable(b"domain", domain);
        self._append_serializable(b"vk", kzg_vk);
    }

    fn set_keyset_commitment(&mut self, keyset_commitment: &SignerSetCommitment) {
        self._append_serializable(b"keyset_commitment", keyset_commitment);
    }

    fn append_public_input(&mut self, public_input: &impl PublicInput) {
        self._append_serializable(b"public_input", public_input);
    }

    fn append_register_commitments(&mut self, register_commitments: &impl RegisterCommitments) {
        self._append_serializable(b"register_commitments", register_commitments);
    }

    fn append_2nd_round_register_commitments(&mut self, register_commitments: &impl RegisterCommitments) {
        self._append_serializable(b"2nd_round_register_commitments", register_commitments);
    }

    fn append_quotient_commitment(&mut self, point: &ark_bw6_761::G1Affine) {
        self._append_serializable(b"quotient", point);
    }

    fn append_register_evaluations(&mut self, evals: &impl RegisterEvaluations) {
        self._append_serializable(b"register_evaluations", evals);
    }

    fn append_quotient_evaluation(&mut self, scalar: &ark_bw6_761::Fr) {
        self._append_serializable(b"quotient_evaluation", scalar);
    }

    fn append_shifted_quotient_evaluation(&mut self, scalar: &ark_bw6_761::Fr) {
        self._append_serializable(b"shifted_quotient_evaluation", scalar);
    }

    fn get_128_bit_challenge(&mut self, label: &'static [u8]) -> ark_bw6_761::Fr;

    fn _append_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize);
}

impl ApkTranscript for Transcript {

    fn get_128_bit_challenge(&mut self, label: &'static [u8]) -> ark_bw6_761::Fr {
        let mut buf = [0u8; 16];
        self.challenge_bytes(label, &mut buf);
        ark_bw6_761::Fr::from_random_bytes(&buf).unwrap()
    }

    fn _append_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize) {
        let mut buffer = vec![0; message.serialized_size()];
        message.serialize(&mut buffer);
        self.append_message(label, &buffer);
    }
}