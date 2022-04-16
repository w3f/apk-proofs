use ark_bls12_377::{Fq as Fr, G1Projective};
use ark_std::rand::Rng;
use ark_std::{test_rng, Zero, One};
use fflonk::pcs::PcsParams;
use merlin::Transcript;

use crate::{Bitmask, Keyset, Proof, Prover, PublicInput, setup, Verifier};
use crate::piop::{RegisterCommitments, RegisterEvaluations};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use ark_std::{end_timer, start_timer, UniformRand};

pub(crate) fn _random_bits<R: Rng>(n: usize, density: f64, rng: &mut R) -> Vec<bool> {
    (0..n).map(|_| rng.gen_bool(density)).collect()
}

pub(crate) fn _random_bitmask<R: Rng>(n: usize, rng: &mut R) -> Vec<Fr> {
    _random_bits(n, 2.0 / 3.0, rng).into_iter()
        .map(|b| if b { Fr::one() } else { Fr::zero() })
        .collect()
}

pub(crate) fn random_pks<R: Rng>(n: usize, rng: &mut R) -> Vec<ark_bls12_377::G1Projective> {
    (0..n)
        .map(|_| G1Projective::rand(rng))
        .collect()
}

fn _test_prove_verify<P, V, PI, E, C, AC>(prove: P, verify: V, log_domain_size: u32, proof_size: usize)
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

    let t_setup = start_timer!(|| "setup");
    // let kzg_params = setup::generate_for_keyset(keyset_size, rng);
    let kzg_params = setup::generate_for_domain(log_domain_size, rng);
    end_timer!(t_setup);

    let pks_commitment_ = start_timer!(|| "signer set commitment");
    let pks_comm = keyset.commit(&kzg_params.ck());
    end_timer!(pks_commitment_);

    let t_prover_new = start_timer!(|| "prover precomputation");
    let prover = Prover::new(
        keyset,
        &pks_comm,
        kzg_params.clone(), //TODO
        Transcript::new(b"apk_proof")
    );
    end_timer!(t_prover_new);

    let verifier = Verifier::new(kzg_params.raw_vk(), pks_comm, Transcript::new(b"apk_proof"));

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


pub fn test_simple_scheme(log_domain_size: u32) {
    _test_prove_verify(
        |prover, bitmask| prover.prove_simple(bitmask),
        |verifier, proof, public_input| verifier.verify_simple(&public_input, &proof),
        log_domain_size,
        (5 * 2 + 6) * 48 // 5C + 6F
    );
}

pub fn test_packed_scheme(log_domain_size: u32) {
    _test_prove_verify(
        |prover, bitmask| prover.prove_packed(bitmask),
        |verifier, proof, public_input| verifier.verify_packed(&public_input, &proof),
        log_domain_size,
        (8 * 2 + 9) * 48 // 8C + 9F
    );
}

pub fn test_counting_scheme(log_domain_size: u32) {
    _test_prove_verify(
        |prover, bitmask| prover.prove_counting(bitmask),
        |verifier, proof, public_input| verifier.verify_counting(&public_input, &proof),
        log_domain_size,
        (7 * 2 + 8) * 48 // 7C + 8F
    );
}