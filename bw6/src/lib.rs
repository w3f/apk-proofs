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

    #[test]
    fn apk_proof() {
        let num_pks = 15;

        let rng = &mut test_rng();

        let signer_set = SignerSet::random(num_pks, rng);

        let setup_ = start_timer!(|| "BW6 setup");
        let params = Params::new(signer_set.size(), rng);
        end_timer!(setup_);

        let pks_domain_size = GeneralEvaluationDomain::<F>::compute_size_of_domain(num_pks).unwrap();

        let pks_commitment_ = start_timer!(|| "signer set commitment");
        let pks_comm = signer_set.commit(&params.get_ck(pks_domain_size));
        end_timer!(pks_commitment_);

        let prover = Prover::new(params.to_pk(), &pks_comm, signer_set.get_all(), Transcript::new(b"apk_proof"));
        let verifier = Verifier::new(params.to_vk(), pks_comm, Transcript::new(b"apk_proof"));

        let b: BitVec = (0..num_pks).map(|_| rng.gen_bool(2.0 / 3.0)).collect();
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