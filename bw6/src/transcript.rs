use merlin::Transcript;
use ark_serialize::CanonicalSerialize;
use ark_ff::{Field, ToBytes};
use ark_ec::ProjectiveCurve;
use crate::signer_set::SignerSetCommitment;
use crate::bls::PublicKey;
use crate::{Bitmask, PublicInput};
use crate::piop::{RegisterCommitments, RegisterEvaluations};

/// E - evaluations
pub(crate) trait ApkTranscript {

    // fn set_protocol_params(&mut self, domain_size: u64, h: &ark_bls12_377::G1Affine);

    fn set_keyset_commitment(&mut self, keyset_commitment: &SignerSetCommitment) {
        self._append_serializable(b"keyset_commitment", keyset_commitment);
    }

    fn append_public_input(&mut self, public_input: & impl PublicInput) {
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
    // fn set_protocol_params(&mut self, domain_size: u64, h: &ark_bls12_377::G1Affine) {
    //     let mut buffer = vec![0; domain_size.serialized_size()];
    //     domain_size.serialize(&mut buffer);
    //     self.append_message(b"domain_size", &buffer);
    //
    //     self._append_bytes(b"h", h);
    // }

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