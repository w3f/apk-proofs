//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

use ark_bw6_761::{BW6_761, Fr};
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::field_new;
use ark_poly::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

pub use bitmask::Bitmask;
pub use keyset::{Keyset, KeysetCommitment};

use crate::kzg::KZG10;
use crate::piop::{RegisterCommitments, RegisterEvaluations};

pub use self::prover::*;
pub use self::verifier::*;
use crate::piop::basic::AffineAdditionEvaluationsWithoutBitmask;
use crate::piop::affine_addition::{PartialSumsCommitments, PartialSumsAndBitmaskCommitments};
use crate::piop::bitmask_packing::{SuccinctAccountableRegisterEvaluations, BitmaskPackingCommitments};
use crate::piop::counting::{CountingEvaluations, CountingCommitments};
use ark_bls12_377::G1Affine;

mod prover;
mod verifier;
pub mod endo;
pub mod utils;

pub mod bls;

mod transcript;

pub mod kzg;
mod fsrng;
mod domains;
mod piop;

pub mod setup;
mod bitmask;
mod keyset;

type UniPoly761 = DensePolynomial<<BW6_761 as PairingEngine>::Fr>;
#[allow(non_camel_case_types)]
type KZG_BW6 = KZG10<BW6_761, UniPoly761>;

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


const H_X: Fr = field_new!(Fr, "0");
const H_Y: Fr = field_new!(Fr, "1");
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

    use ark_std::{end_timer, start_timer, UniformRand};
    use ark_std::convert::TryInto;
    use ark_std::test_rng;
    use ark_std::rand::Rng;
    use merlin::Transcript;
    use ark_ff::{One, Zero};
    use ark_bls12_377::G1Projective;
    use ark_ec::ProjectiveCurve;
    use ark_poly::EvaluationDomain;


    pub fn random_bits<R: Rng>(n: usize, density: f64, rng: &mut R) -> Vec<bool> {
        (0..n).map(|_| rng.gen_bool(density)).collect()
    }

    pub fn random_bitmask<R: Rng>(n: usize, rng: &mut R) -> Vec<Fr> {
        random_bits(n, 2.0 / 3.0, rng).into_iter()
            .map(|b| if b { Fr::one() } else { Fr::zero() })
            .collect()
    }

    pub fn random_pks<R: Rng>(n: usize, rng: &mut R) -> Vec<ark_bls12_377::G1Projective> {
        (0..n)
            .map(|_| G1Projective::rand(rng))
            .collect()
    }


    #[test]
    fn h_is_not_in_g1() {
        let h = point_in_g1_complement();
        assert!(h.is_on_curve());
        assert!(!h.is_in_correct_subgroup_assuming_on_curve());
    }

    fn _test_prove_verify<P, V, PI, E, C, AC>(prove: P, verify: V, proof_size: usize)
        where
            P: Fn(Prover, Bitmask) -> (Proof<E, C, AC>, PI),
            V: Fn(Verifier, Proof<E, C, AC>, PI) -> bool,
            PI: PublicInput,
            E: RegisterEvaluations,
            C: RegisterCommitments,
            AC: RegisterCommitments
    {
        let rng = &mut test_rng();

        let keyset_size = 200;
        let keyset = Keyset::new(random_pks(keyset_size, rng));
        let domain_size = keyset.domain.size();

        let t_setup = start_timer!(|| "setup");
        let kzg_params = setup::generate_for_keyset(keyset_size, rng);
        end_timer!(t_setup);

        let pks_commitment_ = start_timer!(|| "signer set commitment");
        let pks_comm = keyset.commit(&kzg_params.get_pk());
        end_timer!(pks_commitment_);

        let t_prover_new = start_timer!(|| "prover precomputation");
        let prover = Prover::new(
            keyset,
            &pks_comm,
            kzg_params.clone(), //TODO
            Transcript::new(b"apk_proof")
        );
        end_timer!(t_prover_new);

        let verifier = Verifier::new(domain_size, kzg_params.get_vk(), pks_comm, Transcript::new(b"apk_proof"));

        let bits = (0..keyset_size).map(|_| rng.gen_bool(2.0 / 3.0)).collect::<Vec<_>>();
        let b = Bitmask::from_bits(&bits);

        let prove_ = start_timer!(|| "BW6 prove");
        let (proof, public_input) = prove(prover, b.clone());
        end_timer!(prove_);

        let mut serialized_proof = vec![0; proof.serialized_size()];
        proof.serialize(&mut serialized_proof[..]).unwrap();
        let proof = Proof::<E, C, AC>::deserialize(&serialized_proof[..]).unwrap();

        assert_eq!(proof.serialized_size(), proof_size);

        let verify_ = start_timer!(|| "BW6 verify");
        let valid = verify(verifier, proof, public_input);
        end_timer!(verify_);

        assert!(valid);
    }

    #[test]
    fn test_simple_scheme() {
        _test_prove_verify(
            |prover, bitmask| prover.prove_simple(bitmask),
            |verifier, proof, public_input| verifier.verify_simple(&public_input, &proof),
            (5 * 2 + 6) * 48 // 5C + 6F
        );
    }


    #[test]
    fn test_packed_scheme() {
        _test_prove_verify(
            |prover, bitmask| prover.prove_packed(bitmask),
            |verifier, proof, public_input| verifier.verify_packed(&public_input, &proof),
            (8 * 2 + 9) * 48 // 8C + 9F
        );
    }

    #[test]
    fn test_counting_scheme() {
        _test_prove_verify(
            |prover, bitmask| prover.prove_counting(bitmask),
            |verifier, proof, public_input| verifier.verify_counting(&public_input, &proof),
            (7 * 2 + 8) * 48 // 7C + 8F
        );
    }
}