use crate::{PublicKey, KZG_BW6, SecretKey};
use ark_poly_commit::kzg10::Powers;
use ark_bw6_761::{BW6_761, Fr};
use ark_poly::{Evaluations, GeneralEvaluationDomain, EvaluationDomain};
use bitvec::vec::BitVec;
use rand::Rng;
use ark_ec::ProjectiveCurve;
use ark_std::convert::TryInto;
use crate::setup::CommitmentKey;

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

    pub fn commit(&self, ck: CommitmentKey) -> SignerSetCommitment {
        assert!(self.0.len() <= ck.domain.size());

        let (pks_x, pks_y): (Vec<Fr>, Vec<Fr>) = self.0.iter()
            .map(|p| p.0.into_affine())
            .map(|p| (p.x, p.y))
            .unzip();

        let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, ck.domain).interpolate();
        let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, ck.domain).interpolate();

        let pks_x_comm= KZG_BW6::commit(&ck.kzg_ck, &pks_x_poly).unwrap();
        let pks_y_comm= KZG_BW6::commit(&ck.kzg_ck, &pks_y_poly).unwrap();
        SignerSetCommitment {
            pks_x_comm,
            pks_y_comm,
            signer_set_size: self.0.len()
        }
    }

    pub fn random<R: Rng>(num_pks: usize, rng: &mut R) -> Self {
        assert!(num_pks > 1); // https://github.com/arkworks-rs/poly-commit/issues/40
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

    pub fn get_by_mask(&self, b: &BitVec) -> Vec<&PublicKey> {
        self.0.iter().zip(b.iter()).filter(|(_p, b)| **b).map(|(p, _b)| p).collect()
    }
}