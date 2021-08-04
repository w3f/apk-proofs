use apk_proofs::{SignerSet, Setup, SignerSetCommitment, Prover, Verifier, Bitmask, SimpleProof, AccountablePublicInput, hash_to_curve};
use apk_proofs::bls::{PublicKey, SecretKey, Signature};
use apk_proofs::kzg::{VerifierKey, ProverKey};

use ark_serialize::CanonicalSerialize;
use ark_bls12_377::G2Projective;
use ark_bw6_761::BW6_761;
use ark_std::test_rng;

use rand::Rng;
use std::collections::HashSet;
use merlin::Transcript;



// These example sketches the primary intended use case of the crate functionality:
// building communication-efficient light clients for blockchains.

// Here we model a blockchain as a set of validators who are responsible to sign for the chain events.
// The validator set changes in periods of time called 'epochs'. We assume that within an epoch,
// only a fraction of validators in the set is malicious/unresponsive.

// Light client is a resource-constrained blockchain client (think a mobile app or better an Ethereum smart contract),
// that is interested in some of the chain events, but is not able to follow the chain itself.
// Instead it relies on a helper node that provides cryptographic proofs of the events requested by the client
// and doesn't need to be trusted.

// An example of such a proof could be a collection of signatures from the relevant validator set, but it would require
// the client to know all the validators' public keys, that is inefficient. Knowing the aggregate public key
// of the validator set doesn't help, as some of the individual signatures may be missing
// (due to unresponsive/malicious/deluded validators)...



#[derive(Clone)]
struct Validator(SecretKey);

struct Approval {
    comm: SignerSetCommitment,
    sig: Signature,
    pk: PublicKey,
}

impl Validator {
    fn new<R: Rng>(rng: &mut R) -> Self {
        Self(SecretKey::new(rng))
    }

    fn public_key(&self) -> PublicKey {
        (&self.0).into()
    }

    fn approve(&self, new_validator_set: &ValidatorSet, domain_size: usize, kzg_pk: &ProverKey<BW6_761>) -> Approval {
        let new_validator_set_commitment = SignerSet::new(&new_validator_set.public_keys())
            .commit(domain_size, kzg_pk);
        let message = hash_commitment(&new_validator_set_commitment);
        Approval {
            comm: new_validator_set_commitment,
            sig: self.0.sign(&message),
            pk: self.public_key(),
        }
    }
}

#[derive(Clone)]
struct ValidatorSet(Vec<Validator>);

impl ValidatorSet {
    fn new<R: Rng>(size: usize, rng: &mut R) -> Self {
        let validators = (0..size)
            .map(|_| Validator::new(rng))
            .collect();
        Self(validators)
    }

    fn public_keys(&self) -> Vec<PublicKey> {
        self.0.iter()
            .map(|v| v.public_key())
            .collect()
    }

    fn rotate<R: Rng>(&self, domain_size: usize, kzg_pk: &ProverKey<BW6_761>, rng: &mut R) -> (ValidatorSet, Vec<Approval>) {
        let new_validator_set = ValidatorSet::new(self.0.len(), rng);
        let approvals = self.0.iter()
            .filter(|_| rng.gen_bool(0.9))
            .map(|v| v.approve(&new_validator_set, domain_size, kzg_pk))
            .collect();
        (new_validator_set, approvals)
    }
}

struct LightClient {
    domain_size: usize,
    kzg_vk: VerifierKey<BW6_761>,

    current_validator_set_commitment: SignerSetCommitment,
}

impl LightClient {
    fn init(
        domain_size: usize,
        kzg_vk: VerifierKey<BW6_761>,
        genesis_keyset_commitment: SignerSetCommitment,
    ) -> Self {
        Self {
            domain_size,
            kzg_vk,
            current_validator_set_commitment: genesis_keyset_commitment,
        }
    }

    fn verify_aggregates(&mut self,
                         public_input: AccountablePublicInput,
                         proof: &SimpleProof,
                         aggregate_signature: &Signature,
                         new_validator_set_commitment: SignerSetCommitment) {
        let verifier = Verifier::new(self.domain_size, self.kzg_vk.clone(), self.current_validator_set_commitment.clone(), Transcript::new(b"apk_proof"));

        assert!(verifier.verify_simple(&public_input, &proof));
        let aggregate_public_key = public_input.apk;
        let message = hash_commitment(&new_validator_set_commitment);
        assert!(aggregate_public_key.verify(&aggregate_signature, &message));

        self.current_validator_set_commitment = new_validator_set_commitment;
    }
}

struct TrustlessHelper {
    setup: Setup,
    current_validator_set: ValidatorSet,
    prover: Prover
}

impl TrustlessHelper {
    fn new(genesis_validator_set: ValidatorSet, genesis_validator_set_commitment: &SignerSetCommitment, setup: Setup) -> Self {
        let prover = Prover::new(
            &setup,
            genesis_validator_set_commitment,
            genesis_validator_set.public_keys(),
            Transcript::new(b"apk_proof")
        );
        Self {
            setup,
            current_validator_set: genesis_validator_set,
            prover
        }
    }

    fn aggregate_approvals(&mut self, new_validator_set: ValidatorSet, approvals: Vec<Approval>) -> (AccountablePublicInput, SimpleProof, Signature, SignerSetCommitment) {
        let new_validator_set_commitment = &approvals[0].comm;
        let actual_signers = approvals.iter()
            .map(|a| &a.pk)
            .collect::<HashSet<_>>();
        let actual_signers_bitmask = self.current_validator_set.public_keys().iter()
            .map(|pk| actual_signers.contains(pk))
            .collect::<Vec<_>>();

        let (proof, public_input) = self.prover.prove_simple(Bitmask::from_bits(&actual_signers_bitmask));
        let signatures = approvals.iter()
            .map(|a| &a.sig);
        let aggregate_signature = Signature::aggregate(signatures);

        self.current_validator_set = new_validator_set.clone();
        self.prover = Prover::new(
            &self.setup,
            new_validator_set_commitment,
            new_validator_set.public_keys(),
            Transcript::new(b"apk_proof"),
        );

        (public_input, proof, aggregate_signature, new_validator_set_commitment.clone())
    }
}

fn hash_commitment(commitment: &SignerSetCommitment) -> G2Projective {
    let mut buf = vec![0u8; commitment.serialized_size()];
    commitment.serialize(&mut buf[..]).unwrap();
    hash_to_curve(&buf)
}

fn main() {
    let rng = &mut test_rng(); // Don't use in production code!
    let log_keyset_size = 6;
    let keyset_size = 2u64.pow(log_keyset_size) - 1;
    let setup = Setup::generate(log_keyset_size, rng);

    let genesis_validator_set = ValidatorSet::new(keyset_size as usize, rng);
    let genesis_validator_set_commitment = SignerSet::new(&genesis_validator_set.public_keys())
        .commit(setup.domain_size, &setup.kzg_params.get_pk());

    let mut helper = TrustlessHelper::new(genesis_validator_set.clone(), &genesis_validator_set_commitment, setup.clone());
    let mut light_client = LightClient::init(setup.domain_size, setup.kzg_params.get_vk(), genesis_validator_set_commitment);

    let mut current_validator_set = genesis_validator_set;

    for _epoch in 0..2 {
        let (new_validator_set, approvals) = current_validator_set.rotate(setup.domain_size, &setup.kzg_params.get_pk(), rng);

        let (public_input, proof, aggregate_signature, new_validator_set_commitment) =
            helper.aggregate_approvals(new_validator_set.clone(), approvals);

        light_client.verify_aggregates(
            public_input,
            &proof,
            &aggregate_signature,
            new_validator_set_commitment,
        );

        current_validator_set = new_validator_set;
    }
}