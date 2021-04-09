//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

mod prover;
pub use self::prover::*;

mod verifier;
pub use self::verifier::*;

pub mod endo;
pub mod utils;

pub mod bls;

mod transcript;

mod signer_set;
pub use signer_set::{SignerSet, SignerSetCommitment};

mod kzg;
mod fsrng;
mod domains;
mod constraints;
mod piop;

mod setup;
pub use setup::Setup;

mod bitmask;
pub use bitmask::Bitmask;

use ark_poly::univariate::DensePolynomial;
use ark_ec::PairingEngine;

use ark_bw6_761::{BW6_761, Fr, G1Affine};

type UniPoly761 = DensePolynomial<<BW6_761 as PairingEngine>::Fr>;
#[allow(non_camel_case_types)]
type KZG_BW6 = KZG10<BW6_761, UniPoly761>;

use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};

use crate::kzg::KZG10;


#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct BasicRegisterCommitments {
    b_comm: ark_bw6_761::G1Affine,
    acc_comm: (ark_bw6_761::G1Affine, ark_bw6_761::G1Affine),
}


pub trait RegisterCommitments<AC: AdditionalCommitments> {
    fn wrap(
        basic_commitments: BasicRegisterCommitments,
        accountability_commitments: Option<AC>,
    ) -> Self;

    fn get_basic_commitments(&self) -> &BasicRegisterCommitments;

    fn get_additional_commitments(&self) -> Option<AC>;

    fn all_as_vec(&self) -> Vec<ark_bw6_761::G1Affine> {
        let mut res = self.get_basic_commitments().all_as_vec();
        match self.get_additional_commitments() {
            Some(accountable_commitments) => {
                res.extend(accountable_commitments.as_vec());
            },
            _ => {},
        }
        res
    }
}

impl RegisterCommitments<()> for BasicRegisterCommitments {
    fn wrap(basic_commitments: BasicRegisterCommitments, accountability_commitments: Option<()>) -> Self {
        basic_commitments
    }

    fn get_basic_commitments(&self) -> &BasicRegisterCommitments {
        self
    }

    fn get_additional_commitments(&self) -> Option<()> {
        None
    }

    fn all_as_vec(&self) -> Vec<G1Affine> {
        vec![
            self.b_comm,
            self.acc_comm.0,
            self.acc_comm.1,
        ]
    }
}

pub trait AdditionalCommitments {
    fn as_vec(&self) -> Vec<G1Affine>;
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct PackedRegisterCommitments {
    c_comm: ark_bw6_761::G1Affine,
    acc_comm: ark_bw6_761::G1Affine,
}

impl PackedRegisterCommitments {
    pub fn new(c_comm: G1Affine, acc_comm: G1Affine) -> Self {
        PackedRegisterCommitments { c_comm, acc_comm }
    }
}

impl AdditionalCommitments for PackedRegisterCommitments {
    fn as_vec(&self) -> Vec<G1Affine> {
        vec![
            self.c_comm,
            self.acc_comm,
        ]
    }
}

impl AdditionalCommitments for () {
    fn as_vec(&self) -> Vec<G1Affine> {
        unimplemented!()
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct ExtendedRegisterCommitments<AC: CanonicalSerialize + CanonicalDeserialize> {
    basic_commitments: BasicRegisterCommitments,
    additional_commitments: Option<AC>,
}

impl RegisterCommitments<PackedRegisterCommitments> for ExtendedRegisterCommitments<PackedRegisterCommitments> {
    fn wrap(
        basic_commitments: BasicRegisterCommitments,
        packed_commitments: Option<PackedRegisterCommitments>,
    ) -> Self {
        Self {
            basic_commitments,
            additional_commitments: packed_commitments,
        }
    }

    fn get_basic_commitments(&self) -> &BasicRegisterCommitments {
        &self.basic_commitments
    }

    fn get_additional_commitments(&self) -> Option<PackedRegisterCommitments> {
        self.additional_commitments.clone()
    }
}

// #[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E, C> {
    register_commitments: C,
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


use ark_ff::field_new;
const H_X: Fr = field_new!(Fr, "0");
const H_Y: Fr = field_new!(Fr, "1");
fn point_in_g1_complement() -> ark_bls12_377::G1Affine {
    ark_bls12_377::G1Affine::new(H_X, H_Y, false)
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{end_timer, start_timer};
    use merlin::Transcript;
    use ark_std::convert::TryInto;
    use ark_std::test_rng;
    use rand::Rng;
    use crate::constraints::{SuccinctlyAccountableRegisters, Registers, BasicRegisterEvaluations, SuccinctAccountableRegisterEvaluations};

    pub fn random_bits<R: Rng>(size: usize, density: f64, rng: &mut R) -> Vec<bool> {
        (0..size).map(|_| rng.gen_bool(density)).collect()
    }

    #[test]
    fn h_is_not_in_g1() {
        let h = point_in_g1_complement();
        assert!(h.is_on_curve());
        assert!(!h.is_in_correct_subgroup_assuming_on_curve());
    }

    #[test]
    fn test_packed_accountable_scheme() {
        let rng = &mut test_rng();
        let log_domain_size = 10;

        let t_setup = start_timer!(|| "setup");
        let setup = Setup::generate(log_domain_size, rng);
        end_timer!(t_setup);

        let keyset_size = rng.gen_range(0, setup.max_keyset_size()) + 1;
        let keyset_size = keyset_size.try_into().unwrap();
        let signer_set = SignerSet::random(keyset_size, rng);

        let pks_commitment_ = start_timer!(|| "signer set commitment");
        let pks_comm = signer_set.commit(setup.domain_size, &setup.kzg_params.get_pk());
        end_timer!(pks_commitment_);

        let t_prover_new = start_timer!(|| "prover precomputation");
        let prover = Prover::new(
            setup.domain_size,
            setup.kzg_params.get_pk(),
            &pks_comm,
            signer_set.get_all(),
            Transcript::new(b"apk_proof")
        );
        end_timer!(t_prover_new);

        let verifier = Verifier::new(setup.domain_size, setup.kzg_params.get_vk(), pks_comm, Transcript::new(b"apk_proof"));

        let bits = (0..keyset_size).map(|_| rng.gen_bool(2.0 / 3.0)).collect::<Vec<_>>();
        let b = Bitmask::from_bits(&bits);
        let apk = bls::PublicKey::aggregate(signer_set.get_by_mask(&b));

        let prove_ = start_timer!(|| "BW6 prove");
        let proof = prover.prove_packed(&b);
        end_timer!(prove_);

        // let mut serialized_proof = vec![0; proof.serialized_size()];
        // proof.serialize(&mut serialized_proof[..]).unwrap();
        // let proof = Proof::deserialize(&serialized_proof[..]).unwrap();

        let verify_ = start_timer!(|| "BW6 verify");
        let valid = verifier.verify_packed(&apk, &b, &proof);
        end_timer!(verify_);

        assert!(valid);
    }

    #[test]
    fn test_linear_accountable_scheme() {
        let rng = &mut test_rng();
        let log_domain_size = 10;

        let t_setup = start_timer!(|| "setup");
        let setup = Setup::generate(log_domain_size, rng);
        end_timer!(t_setup);

        let keyset_size = rng.gen_range(0, setup.max_keyset_size()) + 1;
        let keyset_size = keyset_size.try_into().unwrap();
        let signer_set = SignerSet::random(keyset_size, rng);

        let pks_commitment_ = start_timer!(|| "signer set commitment");
        let pks_comm = signer_set.commit(setup.domain_size, &setup.kzg_params.get_pk());
        end_timer!(pks_commitment_);

        let t_prover_new = start_timer!(|| "prover precomputation");
        let prover = Prover::new(
            setup.domain_size,
            setup.kzg_params.get_pk(),
            &pks_comm,
            signer_set.get_all(),
            Transcript::new(b"apk_proof")
        );
        end_timer!(t_prover_new);

        let verifier = Verifier::new(setup.domain_size, setup.kzg_params.get_vk(), pks_comm, Transcript::new(b"apk_proof"));

        let bits = (0..keyset_size).map(|_| rng.gen_bool(2.0 / 3.0)).collect::<Vec<_>>();
        let b = Bitmask::from_bits(&bits);
        let apk = bls::PublicKey::aggregate(signer_set.get_by_mask(&b));

        let prove_ = start_timer!(|| "BW6 prove");
        let proof = prover.prove_simple(&b);
        end_timer!(prove_);

        // let mut serialized_proof = vec![0; proof.serialized_size()];
        // proof.serialize(&mut serialized_proof[..]).unwrap();
        // let proof = Proof::deserialize(&serialized_proof[..]).unwrap();

        let verify_ = start_timer!(|| "BW6 verify");
        let valid = verifier.verify_simple(&apk, &b, &proof);
        end_timer!(verify_);

        assert!(valid);
    }
}