use ark_poly::domain::radix2::Radix2EvaluationDomain;
use ark_poly::EvaluationDomain;
use ark_bw6_761::{BW6_761, Fr};
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::{Field, One, batch_inversion};
use rand::Rng;
use crate::kzg::UniversalParams;
use crate::{KZG_BW6, kzg};
use ark_poly_commit::kzg10::Powers;
use ark_std::convert::TryInto;

pub struct Params {
    domain: Radix2EvaluationDomain<Fr>,
    kzg_params: UniversalParams<BW6_761>,

    h: ark_bls12_377::G1Affine,
}

pub struct CommitmentKey<'a> {
    pub domain: Radix2EvaluationDomain<Fr>,
    pub kzg_ck: Powers<'a, BW6_761>,
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
    pub fn generate<R: Rng>(domain_size: u32, rng: &mut R) -> Self {
        assert!(domain_size.is_power_of_two(), "domain size should be a power of 2");
        let n = domain_size.try_into().expect("domain size doesn't fit usize");
        let domain = Radix2EvaluationDomain::<Fr>::new(n).unwrap();
        // deg(q) = 3n-3
        let max_poly_degree = 3 * n - 3; // TODO: assert it fits field's 2-adicity
        let kzg_params = KZG_BW6::setup(max_poly_degree + 1, false, rng).unwrap();

        Self {
            domain,
            kzg_params,
            h: rng.gen::<ark_bls12_377::G1Projective>().into_affine(), //TODO: outside G1
        }
    }

    pub fn max_keyset_size(&self) -> u32 {
        (self.domain.size - 1) as u32
    }

    pub fn get_ck(&self) -> CommitmentKey {
        let n = self.domain.size(); //TODO: smaller cks
        let powers_of_g = self.kzg_params.powers_of_g[..n].to_vec();
        let powers = Powers {
            powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
            powers_of_gamma_g: ark_std::borrow::Cow::default(),
        };
        CommitmentKey {
            domain: self.domain,
            kzg_ck: powers
        }
    }

    pub fn to_pk(&self) -> ProverKey {
        let n = self.domain.size();
        let powers_of_g = self.kzg_params.powers_of_g[..3*n-2].to_vec();
        let powers = Powers {
            powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
            powers_of_gamma_g: ark_std::borrow::Cow::default(),
        };
        ProverKey {
            domain_size: n,
            kzg_ck: powers,
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
        let rng = &mut test_rng();

        let n = 16;
        let params = Params::generate(n, rng);
        let n: usize = n.try_into().unwrap();
        assert_eq!(params.domain.size(), n);

        let z = Fr::rand(rng);
        let evals = params.to_vk().lagrange_evaluations(z);
        assert_eq!(evals.vanishing_polynomial, params.domain.evaluate_vanishing_polynomial(z));
        let coeffs = params.domain.evaluate_all_lagrange_coefficients(z);
        assert_eq!(evals.l_0, coeffs[0]);
        assert_eq!(evals.l_minus_1, coeffs[n - 1]);
    }
}