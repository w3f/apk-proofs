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

mod setup;
pub use setup::Setup;

mod bitmask;
pub use bitmask::Bitmask;

use ark_poly::univariate::DensePolynomial;
use ark_ec::PairingEngine;

use ark_bw6_761::{BW6_761, Fr};

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
    // Prover receives \phi, the constraint polynomials batching challenge, here
    q_comm: ark_bw6_761::G1Affine,
    // Prover receives \zeta, the evaluation point challenge, here
    b_zeta: Fr,
    pks_x_zeta: Fr,
    pks_y_zeta: Fr,
    acc_x_zeta: Fr,
    acc_y_zeta: Fr,
    q_zeta: Fr,
    acc_x_zeta_omega: Fr,
    acc_y_zeta_omega: Fr,
    // Prover receives \nu, the KZG opening batching challenge, here
    w1_proof: ark_bw6_761::G1Affine,
    w2_proof: ark_bw6_761::G1Affine,
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
    use bench_utils::{end_timer, start_timer};
    use merlin::Transcript;
    use ark_std::convert::TryInto;
    use ark_std::test_rng;
    use rand::Rng;

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
    fn apk_proof() {
        let rng = &mut test_rng();
        let log_domain_size = 8;

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