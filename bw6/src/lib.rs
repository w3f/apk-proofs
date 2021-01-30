//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

mod prover;
mod verifier;

pub use self::prover::*;
pub use self::verifier::*;

pub mod endo;
pub mod utils;
pub mod bls;

pub use bls::{Signature, SecretKey, PublicKey};

mod transcript;
mod signer_set;

pub use signer_set::{SignerSet, SignerSetCommitment};

mod kzg;
mod setup;

use setup::{ProverKey, PreparedVerifierKey};

use ark_poly::univariate::DensePolynomial;
use ark_ec::PairingEngine;

use ark_bw6_761::{BW6_761, Fr as F};

type UniPoly761 = DensePolynomial<<BW6_761 as PairingEngine>::Fr>;
#[allow(non_camel_case_types)]
type KZG_BW6 = KZG10<BW6_761, UniPoly761>;


use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};
use crate::kzg::KZG10;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof {
    b_comm: ark_bw6_761::G1Affine,
    acc_x_comm: ark_bw6_761::G1Affine,
    acc_y_comm: ark_bw6_761::G1Affine,
    q_comm: ark_bw6_761::G1Affine,

    w1_proof: ark_bw6_761::G1Affine,
    w2_proof: ark_bw6_761::G1Affine,

    pub b_zeta: F,
    pub pks_x_zeta: F,
    pub pks_y_zeta: F,
    pub acc_x_zeta: F,
    pub acc_y_zeta: F,
    pub acc_x_zeta_omega: F,
    pub acc_y_zeta_omega: F,

    pub q_zeta: F,
}


use ark_ff::field_new;
const H_X: F = field_new!(F, "0");
const H_Y: F = field_new!(F, "1");
fn nums_point_in_g1_complement() -> ark_bls12_377::G1Affine {
    ark_bls12_377::G1Affine::new(H_X, H_Y, false)
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{UniformRand, test_rng};
    use bench_utils::{end_timer, start_timer};
    use merlin::Transcript;
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain};
    use bitvec::vec::BitVec;
    use crate::setup::Params;
    use rand::Rng;
    use ark_std::convert::TryInto;

    #[test]
    fn h_is_not_in_g1() {
        let h = nums_point_in_g1_complement();
        assert!(h.is_on_curve());
        assert!(!h.is_in_correct_subgroup_assuming_on_curve());
    }

    #[test]
    fn apk_proof() {
        let rng = &mut test_rng();

        let log_domain_size = 4;
        let domain_size = 2u32.pow(log_domain_size);

        let t_setup = start_timer!(|| format!("BW6 setup for log(domain_size) = {}", log_domain_size));
        let params = Params::generate(domain_size, rng);
        end_timer!(t_setup);

        let keyset_size = rng.gen_range(1, params.max_keyset_size() + 1);
        let keyset_size = keyset_size.try_into().unwrap();
        let signer_set = SignerSet::random(keyset_size, rng);

        let pks_commitment_ = start_timer!(|| "signer set commitment");
        let pks_comm = signer_set.commit(params.get_ck());
        end_timer!(pks_commitment_);

        let prover = Prover::new(params.to_pk(), &pks_comm, signer_set.get_all(), Transcript::new(b"apk_proof"));
        let verifier = Verifier::new(params.to_vk(), pks_comm, Transcript::new(b"apk_proof"));

        let b: BitVec = (0..keyset_size).map(|_| rng.gen_bool(2.0 / 3.0)).collect();
        let apk = bls::PublicKey::aggregate(signer_set.get_by_mask(&b));

        let prove_ = start_timer!(|| "BW6 prove");
        let proof = prover.prove(&b);
        end_timer!(prove_);

        let mut serialized_proof = vec![0; proof.serialized_size()];
        proof.serialize(&mut serialized_proof[..]).unwrap();
        let proof = Proof::deserialize(&serialized_proof[..]).unwrap();

        let verify_ = start_timer!(|| "BW6 verify");
        let valid = verifier.verify(&apk, &b, &proof);
        end_timer!(verify_);

        assert!(valid);
    }
}