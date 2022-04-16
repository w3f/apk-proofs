//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

use ark_bls12_377::G1Affine;
use ark_bw6_761::{BW6_761, Fr};
use ark_ec::ProjectiveCurve;
use ark_ff::MontFp;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};
use fflonk::pcs::kzg::KZG;

pub use bitmask::Bitmask;
pub use keyset::{Keyset, KeysetCommitment};

use crate::piop::{RegisterCommitments, RegisterEvaluations};
use crate::piop::affine_addition::{PartialSumsAndBitmaskCommitments, PartialSumsCommitments};
use crate::piop::basic::AffineAdditionEvaluationsWithoutBitmask;
use crate::piop::bitmask_packing::{BitmaskPackingCommitments, SuccinctAccountableRegisterEvaluations};
use crate::piop::counting::{CountingCommitments, CountingEvaluations};

pub use self::prover::*;
pub use self::verifier::*;

mod prover;
mod verifier;
pub mod endo;
pub mod utils;

pub mod bls;

mod transcript;

pub mod kzg;
mod fsrng;
pub mod domains;
mod piop;

pub mod setup;
mod bitmask;
mod keyset;
pub mod test_helpers; //TODO: cfgtest

type NewKzgBw6 = KZG<BW6_761>;

// TODO: 1. From trait?
// TODO: 2. remove refs/clones
pub trait PublicInput : CanonicalSerialize + CanonicalDeserialize {
    fn new(apk: &G1Affine, bitmask: &Bitmask) -> Self;
}

// Used in 'basic' and 'packed' schemes
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct AccountablePublicInput {
    pub apk: G1Affine,
    pub bitmask: Bitmask,
}

impl PublicInput for AccountablePublicInput {
    fn new(apk: &G1Affine, bitmask: &Bitmask) -> Self {
        AccountablePublicInput {
            apk: apk.clone(),
            bitmask: bitmask.clone(),
        }
    }
}

// Used in 'counting' scheme
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct CountingPublicInput {
    pub apk: G1Affine,
    pub count: usize,
}

impl PublicInput for CountingPublicInput {
    fn new(apk: &G1Affine, bitmask: &Bitmask) -> Self {
        CountingPublicInput {
            apk: apk.clone(),
            count: bitmask.count_ones(),
        }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: RegisterEvaluations, C: RegisterCommitments, AC: RegisterCommitments> {
    register_commitments: C,
    // 2nd round commitments, used in "packed" scheme after get the bitmask aggregation challenge is received
    additional_commitments: AC,
    // Prover receives \phi, the constraint polynomials batching challenge, here
    q_comm: ark_bw6_761::G1Affine,
    // Prover receives \zeta, the evaluation point challenge, here
    register_evaluations: E,
    q_zeta: Fr,
    r_zeta_omega: Fr,
    // Prover receives \nu, the KZG opening batching challenge, here
    w_at_zeta_proof: ark_bw6_761::G1Affine,
    r_at_zeta_omega_proof: ark_bw6_761::G1Affine,
}

pub type SimpleProof = Proof<AffineAdditionEvaluationsWithoutBitmask, PartialSumsCommitments, ()>;
pub type PackedProof = Proof<SuccinctAccountableRegisterEvaluations, PartialSumsAndBitmaskCommitments, BitmaskPackingCommitments>;
pub type CountingProof = Proof<CountingEvaluations, CountingCommitments, ()>;


const H_X: Fr = MontFp!(Fr, "0");
const H_Y: Fr = MontFp!(Fr, "1");
fn point_in_g1_complement() -> ark_bls12_377::G1Affine {
    ark_bls12_377::G1Affine::new(H_X, H_Y, false)
}

// TODO: switch to better hash to curve when available
pub fn hash_to_curve<G: ProjectiveCurve>(message: &[u8]) -> G {
    use blake2::Digest;
    use ark_std::rand::SeedableRng;

    let seed = blake2::Blake2s::digest(message);
    let rng = &mut rand::rngs::StdRng::from_seed(seed.into());
    G::rand(rng)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers;

    #[test]
    fn h_is_not_in_g1() {
        let h = point_in_g1_complement();
        assert!(h.is_on_curve());
        assert!(!h.is_in_correct_subgroup_assuming_on_curve());
    }

    #[test]
    fn test_simple_scheme() {
        test_helpers::test_simple_scheme(8);
    }


    #[test]
    fn test_packed_scheme() {
        test_helpers::test_packed_scheme(8);
    }

    #[test]
    fn test_counting_scheme() {
        test_helpers::test_counting_scheme(8);
    }
}