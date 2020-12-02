//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

mod prover;
mod verifier;

pub use self::prover::*;
pub use self::verifier::*;

pub mod endo;
pub mod utils;
pub mod bls;
pub use bls::{Signature, SecretKey, PublicKey};

use ark_ff::{One, Field, batch_inversion};
use ark_poly::{Evaluations, EvaluationDomain, GeneralEvaluationDomain, Radix2EvaluationDomain};
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::kzg10::{KZG10, Powers};
use ark_ec::{ProjectiveCurve, PairingEngine};

use bitvec::vec::BitVec;
use rand::Rng;
use ark_bw6_761::{BW6_761, Fr as F};

type UniPoly_761 = DensePolynomial<<BW6_761 as PairingEngine>::Fr>;
type KZG_BW6 = KZG10<BW6_761, UniPoly_761>;


pub struct Params {
    domain: Radix2EvaluationDomain<F>,
    kzg_params: ark_poly_commit::kzg10::UniversalParams<BW6_761>,

    h: ark_bls12_377::G1Affine,
}

pub struct ProverKey<'a> {
    domain_size: usize,
    kzg_ck: Powers<'a, BW6_761>,
    domain: Radix2EvaluationDomain<F>,
    h: ark_bls12_377::G1Affine,
}

pub struct VerifierKey {
    domain_size: u64,
    kzg_vk: ark_poly_commit::kzg10::VerifierKey<BW6_761>,
    domain: Radix2EvaluationDomain<F>,
    h: ark_bls12_377::G1Affine,
}

pub struct LagrangeEvaluations {
    vanishing_polynomial: F,
    l_0: F,
    l_minus_1: F,
}

impl VerifierKey {
    pub fn lagrange_evaluations(&self, zeta: F) -> LagrangeEvaluations {
        let mut zeta_n = zeta;
        for _ in 0..self.domain.log_size_of_group {
            zeta_n.square_in_place();
        }
        assert_eq!(zeta_n, zeta.pow([self.domain_size]));
        let zeta_n_minus_one= zeta_n - F::one();
        let zeta_n_minus_one_div_n = zeta_n_minus_one * self.domain.size_inv;

        let mut inv = [zeta - F::one(), self.domain.group_gen * zeta - F::one()];
        batch_inversion(&mut inv);
        LagrangeEvaluations {
            vanishing_polynomial: zeta_n_minus_one,
            l_0: zeta_n_minus_one_div_n * inv[0],
            l_minus_1: zeta_n_minus_one_div_n * inv[1],
        }
    }
}

impl Params {
    pub fn new<R: Rng>(max_pks: usize, rng: &mut R) -> Self {
        let min_domain_size = max_pks + 1; // to initialize the acc with h
        let domain = Radix2EvaluationDomain::<F>::new(min_domain_size).unwrap();
        let n = domain.size();

        // deg(q) = 3n-3
        let kzg_params = KZG_BW6::setup(3*n-2, false, rng).unwrap();

        Self {
            domain,
            kzg_params,

            h: rng.gen::<ark_bls12_377::G1Projective>().into_affine(), //TODO: outside G1
        }
    }

    pub fn get_ck(&self, m: usize) -> Powers<BW6_761> {
        let powers_of_g = self.kzg_params.powers_of_g[..m].to_vec();
        Powers {
            powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
            powers_of_gamma_g: ark_std::borrow::Cow::default(),
        }
    }

    pub fn to_pk(&self) -> ProverKey {
        let n = self.domain.size();
        ProverKey {
            domain_size: n,
            kzg_ck: self.get_ck(3*n-2),
            domain: self.domain,
            h: self.h
        }
    }

    pub fn to_vk(&self) -> VerifierKey {
        let vk = ark_poly_commit::kzg10::VerifierKey {
            g: self.kzg_params.powers_of_g[0],
            gamma_g: self.kzg_params.powers_of_gamma_g[&0],
            h: self.kzg_params.h,
            beta_h: self.kzg_params.beta_h,
            prepared_h: self.kzg_params.prepared_h.clone(),
            prepared_beta_h: self.kzg_params.prepared_beta_h.clone(),
        };
        VerifierKey {
            domain_size: self.domain.size,
            kzg_vk: vk,
            domain: self.domain,
            h: self.h
        }
    }
}

struct SignerSet(Vec<PublicKey>);

impl SignerSet {
    pub fn size(&self) -> usize {
        self.0.len()
    }

    pub fn commit(&self, ck: &Powers<BW6_761>) -> (ark_bw6_761::G1Affine, ark_bw6_761::G1Affine) {
        let m = self.0.len();
        // assert_eq!(m, powers.len());
        // as now we use ifft to compute the polynomials, we require
        let domain = GeneralEvaluationDomain::<F>::new(m).unwrap();
        assert_eq!(domain.size(), ck.powers_of_g.len());

        let (pks_x, pks_y): (Vec<F>, Vec<F>) = self.0.iter()
            .map(|p| p.0.into_affine())
            .map(|p| (p.x, p.y))
            .unzip();

        let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, domain).interpolate();
        let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, domain).interpolate();

        let (pks_x_comm, _) = KZG_BW6::commit(ck, &pks_x_poly, None, None).unwrap();
        let (pks_y_comm, _) = KZG_BW6::commit(ck, &pks_y_poly, None, None).unwrap();
        (pks_x_comm.0, pks_y_comm.0)
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

    pub zeta: F,
    pub phi: F,
    pub nu: F,
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;
    use ark_ff::{test_rng, UniformRand};

    #[test]
    fn test_pks_commitment() {
        let num_pks = 10;

        let rng = &mut test_rng();

        let signer_set = SignerSet::random(num_pks, rng);

        let params = Params::new(signer_set.size(), rng);

        let pks_domain_size = GeneralEvaluationDomain::<F>::compute_size_of_domain(num_pks).unwrap();
        let (pks_x_comm, pks_y_comm) = signer_set.commit(&params.get_ck(pks_domain_size));

        let b: BitVec = (0..num_pks).map(|_| rng.gen_bool(2.0 / 3.0)).collect();

        let apk = bls::PublicKey::aggregate(signer_set.get_by_mask(&b));

        let proof = prove(&b, signer_set.get_all(), &params.to_pk());
        assert!(verify(&pks_x_comm, &pks_y_comm, apk, &b, &proof, &params.to_vk()));
    }

    #[test]
    fn test_lagrange_evaluations() {
        let n = 16;
        let rng = &mut test_rng();
        let params = Params::new(n-1, rng);
        assert_eq!(params.domain.size(), n);

        let z = F::rand(rng);
        let evals = params.to_vk().lagrange_evaluations(z);
        assert_eq!(evals.vanishing_polynomial, params.domain.evaluate_vanishing_polynomial(z));
        let coeffs =  params.domain.evaluate_all_lagrange_coefficients(z);
        assert_eq!(evals.l_0, coeffs[0]);
        assert_eq!(evals.l_minus_1, coeffs[n-1]);
    }
}