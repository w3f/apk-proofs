//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

use ark_bls12_377::{
    g1::Config as G1Config, g2::Config as G2Config, G1Affine, G1Projective, G2Projective,
};
use ark_bw6_761::{Fr, BW6_761};
use ark_ec::{
    hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher, HashToCurve},
    models::short_weierstrass::Projective,
};
use ark_ff::{field_hashers::DefaultFieldHasher, MontFp};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use fflonk::pcs::kzg::KZG;
use sha2::Sha256;

pub use bitmask::Bitmask;
pub use keyset::{Keyset, KeysetCommitment};

use crate::piop::affine_addition::{PartialSumsAndBitmaskCommitments, PartialSumsCommitments};
use crate::piop::basic::AffineAdditionEvaluationsWithoutBitmask;
use crate::piop::bitmask_packing::{
    BitmaskPackingCommitments, SuccinctAccountableRegisterEvaluations,
};
use crate::piop::counting::{CountingCommitments, CountingEvaluations};
use crate::piop::{RegisterCommitments, RegisterEvaluations};

pub use self::prover::*;
pub use self::verifier::*;

pub mod endo;
mod prover;
pub mod utils;
mod verifier;

pub mod bls;

mod transcript;

pub mod domains;
mod fsrng;
mod piop;

mod bitmask;
mod keyset;
pub mod setup;
pub mod test_helpers; //TODO: cfgtest

type NewKzgBw6 = KZG<BW6_761>;

pub const DST_G1: &[u8] = b"APK-PROOF-with-BLS12377G1_XMD:SHA-256_SSWU_RO_";
pub const DST_G2: &[u8] = b"APK-PROOF-with-BLS12377G2_XMD:SHA-256_SSWU_RO_";

// TODO: 1. From trait?
// TODO: 2. remove refs/clones
pub trait PublicInput: CanonicalSerialize + CanonicalDeserialize {
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
pub type PackedProof = Proof<
    SuccinctAccountableRegisterEvaluations,
    PartialSumsAndBitmaskCommitments,
    BitmaskPackingCommitments,
>;
pub type CountingProof = Proof<CountingEvaluations, CountingCommitments, ()>;

const H_X: Fr = MontFp!("0");
const H_Y: Fr = MontFp!("1");
fn point_in_g1_complement() -> ark_bls12_377::G1Affine {
    ark_bls12_377::G1Affine::new_unchecked(H_X, H_Y)
}

pub fn hash_to_curve_g1(message: &[u8]) -> G1Projective {
    let wb_to_curve_hasher = MapToCurveBasedHasher::<
        Projective<G1Config>,
        DefaultFieldHasher<Sha256>,
        WBMap<G1Config>,
    >::new(DST_G1)
    .unwrap();
    wb_to_curve_hasher.hash(message).unwrap().into()
}

pub fn hash_to_curve_g2(message: &[u8]) -> G2Projective {
    let wb_to_curve_hasher = MapToCurveBasedHasher::<
        Projective<G2Config>,
        DefaultFieldHasher<Sha256>,
        WBMap<G2Config>,
    >::new(DST_G2)
    .unwrap();
    wb_to_curve_hasher.hash(message).unwrap().into()
}

#[cfg(test)]
mod tests {
    use crate::test_helpers;

    use super::*;

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
