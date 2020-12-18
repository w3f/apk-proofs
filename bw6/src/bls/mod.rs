use rand::Rng;
use std::borrow::Borrow;
use std::ops::Neg;

use ark_ff::{One, UniformRand};
use ark_ec::{AffineCurve, ProjectiveCurve, PairingEngine};
use ark_bls12_377::{G2Projective, Fr, G1Projective, Bls12_377, G1Affine, Fq12};
use ark_serialize::*;



#[derive(Clone, Debug)]
pub struct Signature(G2Projective);

impl From<G2Projective> for Signature {
    fn from(sig: G2Projective) -> Signature {
        Signature(sig)
    }
}

impl AsRef<G2Projective> for Signature {
    fn as_ref(&self) -> &G2Projective {
        &self.0
    }
}

impl Signature {
    pub fn aggregate<S: Borrow<Signature>>(signatures: impl IntoIterator<Item=S>) -> Signature {
        signatures
            .into_iter()
            .map(|s| s.borrow().0)
            .sum::<G2Projective>()
            .into()
    }
}



#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct SecretKey(Fr);

impl From<Fr> for SecretKey {
    fn from(sk: Fr) -> SecretKey {
        SecretKey(sk)
    }
}

impl AsRef<Fr> for SecretKey {
    fn as_ref(&self) -> &Fr {
        &self.0
    }
}

impl SecretKey {
    pub fn new<R: Rng>(rng: &mut R) -> SecretKey {
        SecretKey(Fr::rand(rng))
    }

    pub fn sign(&self, message: &G2Projective) -> Signature {
	let mut message = message.into_affine(); 
        message.mul(*self.as_ref()).into()
    }
}



#[derive(Clone, Debug)]
pub struct PublicKey(pub G1Projective); // TODO: remove pub

impl From<G1Projective> for PublicKey {
    fn from(pk: G1Projective) -> PublicKey {
        PublicKey(pk)
    }
}

impl From<&SecretKey> for PublicKey {
    fn from(sk: &SecretKey) -> PublicKey {
	
        G1Projective::prime_subgroup_generator().into_affine().mul(*sk.as_ref()).into()
    }
}

impl PublicKey {
    pub fn aggregate<P: Borrow<PublicKey>>(public_keys: impl IntoIterator<Item = P>) -> PublicKey {
        public_keys
            .into_iter()
            .map(|s| s.borrow().0)
            .sum::<G1Projective>()
            .into()
    }

    pub fn verify(&self, signature: &Signature, message: &G2Projective) -> bool {
        Bls12_377::product_of_pairings(&vec![
            (
                G1Affine::prime_subgroup_generator().neg().into(),
                signature.as_ref().into_affine().into(),
            ),
            (
                self.0.into_affine().into(),
                message.into_affine().into(),
            ),
        ]) == Fq12::one()
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::test_rng;

    #[test]
    fn test_apk() {
        let rng = &mut test_rng();
        let message = G2Projective::rand(rng);

        let sks = (0..10).map(|_| SecretKey::new(rng)).collect::<Vec<_>>();
        let pks = sks.iter().map(PublicKey::from).collect::<Vec<_>>();
        let sigs = sks.iter().map(|sk| sk.sign(&message)).collect::<Vec<_>>();
        pks.iter().zip(sigs.iter()).for_each(|(pk, sig)| assert!(pk.verify(sig, &message)));

        let apk = PublicKey::aggregate(pks);
        let asig = Signature::aggregate(sigs);
        assert!(apk.verify(&asig, &message));
    }
}
