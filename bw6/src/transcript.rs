use merlin::Transcript;
use ark_serialize::CanonicalSerialize;
use ark_ff::{Field, ToBytes};
use crate::PublicKey;
use ark_ec::ProjectiveCurve;
use bitvec::vec::BitVec;
use crate::signer_set::SignerSetCommitment;

pub trait ApkTranscript {

    fn set_protocol_params(&mut self, domain_size: u64, h: &ark_bls12_377::G1Affine);

    fn set_signer_set(&mut self, signer_set_comm: &SignerSetCommitment) {
        self._append_bytes(b"pks_x_comm", &signer_set_comm.pks_x_comm);
        self._append_bytes(b"pks_y_comm", &signer_set_comm.pks_y_comm);
        self._append_bytes(b"pks_size", &signer_set_comm.signer_set_size);
    }

    fn append_public_input(&mut self, apk: &PublicKey, bitmask: &BitVec);

    fn append_proof_point(&mut self, label: &'static [u8], point: &ark_bw6_761::G1Affine) {
        self._append_bytes(label, point);
    }

    fn append_proof_scalar(&mut self, label: &'static [u8], scalar: &ark_bw6_761::Fr) {
        self._append_bytes(label, scalar);
    }

    fn get_128_bit_challenge(&mut self, label: &'static [u8]) -> ark_bw6_761::Fr;

    fn _append_bytes<T: ToBytes>(&mut self, label: &'static [u8], message: &T);
}

impl ApkTranscript for Transcript {
    fn set_protocol_params(&mut self, domain_size: u64, h: &ark_bls12_377::G1Affine) {
        let mut buffer = vec![0; domain_size.serialized_size()];
        domain_size.serialize(&mut buffer);
        self.append_message(b"domain_size", &buffer);

        self._append_bytes(b"h", h);
    }

    fn append_public_input(&mut self, apk: &PublicKey, bitmask: &BitVec) {
        let apk = apk.0.into_affine();
        self._append_bytes(b"apk", &apk);

        let bitmask = bitmask.clone().into_vec();
        let mut buffer = vec![0; bitmask.serialized_size()];
        bitmask.serialize(&mut buffer);
        self.append_message(b"bitmask", &buffer);
    }

    fn get_128_bit_challenge(&mut self, label: &'static [u8]) -> ark_bw6_761::Fr {
        let mut buf = [0u8; 16];
        self.challenge_bytes(label, &mut buf);
        ark_bw6_761::Fr::from_random_bytes(&buf).unwrap()
    }

    fn _append_bytes<T: ToBytes>(&mut self, label: &'static [u8], message: &T) {
        let mut buf = Vec::new(); //TODO: suboptimal
        message.write(&mut buf);
        self.append_message(label, &buf);
    }
}