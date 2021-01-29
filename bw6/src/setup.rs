use ark_poly::domain::radix2::Radix2EvaluationDomain;
use ark_poly::EvaluationDomain;
use ark_bw6_761::{BW6_761, Fr};
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::{Field, One, batch_inversion};
use rand::Rng;
use crate::kzg::UniversalParams;
use crate::{KZG_BW6, kzg};
use ark_poly_commit::kzg10::Powers;

pub struct Params {
    domain: Radix2EvaluationDomain<Fr>,
    kzg_params: UniversalParams<BW6_761>,

    h: ark_bls12_377::G1Affine,
}

pub struct ProverKey<'a> {
    pub domain_size: usize,
    pub kzg_ck: Powers<'a, BW6_761>,
    pub domain: Radix2EvaluationDomain<Fr>,
    pub h: ark_bls12_377::G1Affine,
}

pub struct PreparedVerifierKey {
    pub domain_size: u64,
    pub domain: Radix2EvaluationDomain<Fr>,
    pub h: ark_bls12_377::G1Affine,

    pub kzg_vk_prepared: kzg::PreparedVerifierKey<BW6_761>,
}

pub struct LagrangeEvaluations {
    pub vanishing_polynomial: Fr,
    pub l_0: Fr,
    pub l_minus_1: Fr,
}

impl PreparedVerifierKey {
    pub fn lagrange_evaluations(&self, zeta: Fr) -> LagrangeEvaluations {
        let mut zeta_n = zeta;
        for _ in 0..self.domain.log_size_of_group {
            zeta_n.square_in_place();
        }
        assert_eq!(zeta_n, zeta.pow([self.domain_size]));
        let zeta_n_minus_one = zeta_n - Fr::one();
        let zeta_n_minus_one_div_n = zeta_n_minus_one * self.domain.size_inv;

        let mut inv = [zeta - Fr::one(), self.domain.group_gen * zeta - Fr::one()];
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
        let domain = Radix2EvaluationDomain::<Fr>::new(min_domain_size).unwrap();
        let n = domain.size();

        // deg(q) = 3n-3
        let kzg_params = KZG_BW6::setup(3 * n - 2, false, rng).unwrap();

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
            kzg_ck: self.get_ck(3 * n - 2),
            domain: self.domain,
            h: self.h,
        }
    }

    pub fn to_vk(&self) -> PreparedVerifierKey {
        PreparedVerifierKey {
            domain_size: self.domain.size,
            domain: self.domain,
            h: self.h,

            kzg_vk_prepared: kzg::PreparedVerifierKey {
                g: self.kzg_params.powers_of_g[0],
                prepared_h: self.kzg_params.h.into(),
                prepared_beta_h: self.kzg_params.beta_h.into(),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{UniformRand, test_rng};

    #[test]
    fn test_lagrange_evaluations() {
        let n = 16;
        let rng = &mut test_rng();
        let params = Params::new(n - 1, rng);
        assert_eq!(params.domain.size(), n);

        let z = Fr::rand(rng);
        let evals = params.to_vk().lagrange_evaluations(z);
        assert_eq!(evals.vanishing_polynomial, params.domain.evaluate_vanishing_polynomial(z));
        let coeffs = params.domain.evaluate_all_lagrange_coefficients(z);
        assert_eq!(evals.l_0, coeffs[0]);
        assert_eq!(evals.l_minus_1, coeffs[n - 1]);
    }
}