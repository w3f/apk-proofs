use merlin::{Transcript, TranscriptRng};
use rand::{CryptoRng, Error, RngCore};

struct DummyRng;

impl RngCore for DummyRng {
    fn next_u32(&mut self) -> u32 {
        unimplemented!()
    }

    fn next_u64(&mut self) -> u64 {
        unimplemented!()
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        dest.iter_mut().for_each(|byte| *byte = 0u8);
    }

    fn try_fill_bytes(&mut self, _dest: &mut [u8]) -> Result<(), Error> {
        unimplemented!()
    }
}

impl CryptoRng for DummyRng {}

pub fn fiat_shamir_rng(transcript: &mut Transcript) -> TranscriptRng {
    transcript
        .build_rng()
        .rekey_with_witness_bytes(b"verifier_secret", &[42])//TODO: Does verifier know secrets?
        .finalize(&mut DummyRng)
}