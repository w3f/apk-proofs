use crate::{hash_to_curve, KZG_BW6};
use ark_poly::{Radix2EvaluationDomain, EvaluationDomain, Evaluations};
use ark_bls12_377::{G1Affine, G1Projective};
use ark_bw6_761::{Fr, BW6_761};
use ark_poly::univariate::DensePolynomial;
use ark_ec::ProjectiveCurve;
use crate::kzg::ProverKey;
use crate::domains::Domains;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct KeysetCommitment {
    pub pks_comm: (ark_bw6_761::G1Affine, ark_bw6_761::G1Affine),
    pub keyset_size: usize
}

#[derive(Clone)]
pub struct Keyset {
    pub pks: Vec<G1Projective>,
    pub domain: Radix2EvaluationDomain<Fr>,
    pub pks_polys: [DensePolynomial<Fr>; 2],
    pub pks_evals_x4: Option<[Evaluations<Fr, Radix2EvaluationDomain<Fr>>; 2]>,
}

impl Keyset {
    pub fn new(pks: Vec<G1Projective>) -> Self {
        let min_domain_size = pks.len() + 1; // extra 1 accounts apk accumulator initial value
        let domain = Radix2EvaluationDomain::<Fr>::new(min_domain_size).unwrap();
        let mut padded_pks = pks.clone();
        padded_pks.resize_with(domain.size(), || hash_to_curve::<ark_bls12_377::G1Projective>(b"apk-proofs"));

        // convert into affine coordinates to commit
        let (pks_x, pks_y) = G1Projective::batch_normalization_into_affine(&padded_pks).iter()
            .map(|p| (p.x, p.y))
            .unzip();
        let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, domain).interpolate();
        let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, domain).interpolate();
        Self {
            pks: pks,
            domain,
            pks_polys: [pks_x_poly, pks_y_poly],
            pks_evals_x4: None,
        }
    }

    // actual number of signers, not including the padding
    pub fn size(&self) -> usize {
        self.pks.len()
    }

    pub fn amplify(&mut self) {
        let domains = Domains::new(self.domain.size());
        let pks_evals_x4 = self.pks_polys.clone().map(|z| domains.amplify_polynomial(&z));
        self.pks_evals_x4 = Some(pks_evals_x4);
    }

    pub fn commit(&self, kzg_pk: &ProverKey<BW6_761>) -> KeysetCommitment {
        assert!(self.domain.size() <= kzg_pk.max_coeffs());
        let pks_x_comm= KZG_BW6::commit(kzg_pk, &self.pks_polys[0]);
        let pks_y_comm= KZG_BW6::commit(kzg_pk, &self.pks_polys[1]);
        KeysetCommitment {
            pks_comm: (pks_x_comm, pks_y_comm),
            keyset_size: self.pks.len(),
        }
    }

    pub fn aggregate(&self, bitmask: &[bool]) -> ark_bls12_377::G1Projective {
        assert_eq!(bitmask.len(), self.size());
        bitmask.iter()
            .zip(self.pks.iter())
            .filter(|(b, _p)| **b)
            .map(|(_b, p)| p)
            .sum()
    }
}
