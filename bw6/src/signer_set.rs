use crate::{KZG_BW6, Bitmask};
use ark_bw6_761::{Fr, BW6_761};
use ark_poly::{Evaluations, EvaluationDomain, Radix2EvaluationDomain};
use ark_ec::ProjectiveCurve;
use rand::Rng;


use crate::kzg::ProverKey;
use crate::bls::{PublicKey, SecretKey};

pub struct SignerSet(Vec<PublicKey>);

pub struct SignerSetCommitment {
    pub pks_x_comm: ark_bw6_761::G1Affine,
    pub pks_y_comm: ark_bw6_761::G1Affine,
    pub signer_set_size: usize
}

impl SignerSet {
    pub fn size(&self) -> usize {
        self.0.len()
    }

    pub fn commit(&self, domain_size: usize, kzg_pk: &ProverKey<BW6_761>) -> SignerSetCommitment {
        assert!(domain_size.is_power_of_two());
        assert!(self.size() + 1 <= domain_size); // accounts for accumulator initial value h
        assert!(domain_size <= kzg_pk.max_coeffs());

        let (pks_x, pks_y): (Vec<Fr>, Vec<Fr>) = self.0.iter()
            .map(|p| p.0.into_affine())
            .map(|p| (p.x, p.y))
            .unzip();

        let domain = Radix2EvaluationDomain::<Fr>::new(domain_size).unwrap();
        let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, domain).interpolate();
        let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, domain).interpolate();

        let pks_x_comm= KZG_BW6::commit(kzg_pk, &pks_x_poly);
        let pks_y_comm= KZG_BW6::commit(kzg_pk, &pks_y_poly);
        SignerSetCommitment {
            pks_x_comm,
            pks_y_comm,
            signer_set_size: self.0.len()
        }
    }

    pub fn random<R: Rng>(num_pks: usize, rng: &mut R) -> Self {
        assert!(num_pks > 0);
        Self(
            (0..num_pks)
                .map(|_| SecretKey::new(rng))
                .map(|sk| PublicKey::from(&sk))
                .collect::<Vec<_>>()
        )
    }

    pub fn get_all(&self) -> &[PublicKey] {
        return self.0.as_slice();
    }

    pub fn get_by_mask(&self, b: &Bitmask) -> Vec<&PublicKey> {
        self.0.iter().zip(b.to_bits().iter()).filter(|(_p, b)| **b).map(|(p, _b)| p).collect()
    }
}