use ark_bw6_761::{Fr, G1Affine};
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

use crate::{Bitmask, utils};
use crate::domains::Domains;
use crate::utils::LagrangeEvaluations;

pub mod packed;
pub mod affine_addition;
pub mod basic;
pub mod bit_packing;

pub trait RegisterCommitments: CanonicalSerialize + CanonicalDeserialize {
    fn as_vec(&self) -> Vec<G1Affine>;
}

pub trait RegisterPolynomials {
    type C: RegisterCommitments;
    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> Self::C;
}

impl RegisterCommitments for () {
    fn as_vec(&self) -> Vec<G1Affine> {
        vec![]
    }
}

impl RegisterPolynomials for () {
    type C = ();

    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, _f: F) -> Self::C {
        ()
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct PackedRegisterCommitments {
    pub c_comm: ark_bw6_761::G1Affine,
    pub acc_comm: ark_bw6_761::G1Affine,
}

impl PackedRegisterCommitments {
    pub fn new(c_comm: G1Affine, acc_comm: G1Affine) -> Self {
        PackedRegisterCommitments { c_comm, acc_comm }
    }
}

impl RegisterCommitments for PackedRegisterCommitments {
    fn as_vec(&self) -> Vec<G1Affine> {
        vec![
            self.c_comm,
            self.acc_comm,
        ]
    }
}

pub trait Protocol {
    type P1: RegisterPolynomials;
    type P2: RegisterPolynomials;
    type E: RegisterEvaluations;

    fn init(domains: Domains, bitmask: Bitmask, pks: Vec<ark_bls12_377::G1Affine>) -> Self;

    fn get_1st_round_register_polynomials(&self) -> Self::P1;

    fn get_2nd_round_register_polynomials(&mut self, verifier_challenge: Fr) -> Self::P2;


    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>>;
    //TODO: remove domains param
    fn compute_quotient_polynomial(&self, phi: Fr, domains: &Domains) -> DensePolynomial<Fr> {
        let w = utils::randomize(phi, &self.compute_constraint_polynomials());
        let (q_poly, r) = domains.compute_quotient(&w);
        assert_eq!(r, DensePolynomial::zero());
        q_poly
    }

    // TODO: move zeta_minus_omega_inv param to evaluations
    fn evaluate_register_polynomials(&self, point: Fr) -> Self::E;
    // TODO: move zeta_minus_omega_inv param to evaluations
    fn compute_linearization_polynomial(&self, evaluations: &Self::E, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr>;

    fn get_all_register_polynomials(self) -> Vec<DensePolynomial<Fr>>;
}

pub struct BitmaskPackingPolynomials {
    pub c_poly: DensePolynomial<Fr>,
    pub acc_poly: DensePolynomial<Fr>,
}

impl BitmaskPackingPolynomials {
    pub fn new(c_poly: DensePolynomial<Fr>, acc_poly: DensePolynomial<Fr>) -> Self {
        BitmaskPackingPolynomials { c_poly, acc_poly }
    }

    //TODO: &self
    pub fn to_vec(self) -> Vec<DensePolynomial<Fr>> {
        vec![
            self.c_poly,
            self.acc_poly,
        ]
    }
}

impl RegisterPolynomials for BitmaskPackingPolynomials {
    type C = PackedRegisterCommitments;

    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> Self::C {
            PackedRegisterCommitments::new(
                f(&self.c_poly),
                f(&self.acc_poly),
            )
    }
}

pub trait RegisterEvaluations: CanonicalSerialize + CanonicalDeserialize {
    type AC: RegisterCommitments;
    type C: RegisterCommitments;

    fn as_vec(&self) -> Vec<Fr>;
    fn get_bitmask(&self) -> Fr;
    fn restore_commitment_to_linearization_polynomial(
        &self,
        phi: Fr,
        zeta_minus_omega_inv: Fr,
        commitments: &Self::C,
        extra_commitments: &Self::AC,
    ) -> ark_bw6_761::G1Projective;

    fn evaluate_constraint_polynomials(
        &self,
        apk: ark_bls12_377::G1Affine,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
        r: Fr,
        bitmask: &Bitmask,
        domain_size: u64,
    ) -> Vec<Fr>;

    //TODO: move somewhere
    fn is_accountable(&self) -> bool;
}
