use ark_poly::domain::radix2::Radix2EvaluationDomain;
use ark_poly::EvaluationDomain;
use ark_bw6_761::{BW6_761, Fr};
use ark_ec::ProjectiveCurve;
use ark_ff::{Field, One, batch_inversion};
use rand::Rng;
use crate::{KZG_BW6, kzg};
use ark_poly_commit::kzg10::Powers;
use ark_std::convert::TryInto;
use crate::kzg::VerifierKey;

pub struct Params {
    domain: Radix2EvaluationDomain<Fr>,
    kzg_params: kzg::Params<BW6_761>,

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

impl Params {
    pub fn generate<R: Rng>(domain_size: u64, rng: &mut R) -> Self {
        assert!(domain_size.is_power_of_two(), "domain size should be a power of 2");
        let n = domain_size.try_into().expect("domain size doesn't fit usize");
        let domain = Radix2EvaluationDomain::<Fr>::new(n).unwrap();
        // deg(q) = 3n-3
        let max_poly_degree = 3 * n - 3; // TODO: assert it fits field's 2-adicity
        let kzg_params = KZG_BW6::setup(max_poly_degree + 1, rng).unwrap();

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

    pub fn to_vk(&self) -> VerifierKey<BW6_761> {
        self.kzg_params.get_vk()
    }
}