//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

pub mod bls;
pub use bls::{Signature, SecretKey, PublicKey};