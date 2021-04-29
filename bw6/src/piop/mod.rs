use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_bw6_761::{Fr, G1Affine};

use crate::domains::Domains;
use crate::{utils, Bitmask};
use crate::constraints::Registers;
use ark_bls12_377::Fq;

use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};

pub mod packed;

pub trait RegisterCommitments {
    fn as_vec(&self) -> Vec<G1Affine>;
}

pub trait RegisterPolys {
    type C: RegisterCommitments;
    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> Self::C;
}

impl RegisterCommitments for () {
    fn as_vec(&self) -> Vec<G1Affine> {
        vec![]
    }
}

impl RegisterPolys for () {
    type C = ();

    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> Self::C {
        ()
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct PartialSumsCommitments(pub ark_bw6_761::G1Affine, pub ark_bw6_761::G1Affine);

impl RegisterCommitments for PartialSumsCommitments {
    fn as_vec(&self) -> Vec<G1Affine> {
        vec![
            self.0,
            self.1,
        ]
    }
}

pub struct PartialSumsPolynomials(pub DensePolynomial<Fr>, pub DensePolynomial<Fr>);

impl RegisterPolys for PartialSumsPolynomials {
    type C = PartialSumsCommitments;

    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> PartialSumsCommitments {
        PartialSumsCommitments(f(&self.0), f(&self.1))
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

pub trait Protocol<E> {
    type P1: RegisterPolys;
    type P2: RegisterPolys;

    fn init(domains: Domains, bitmask: &Bitmask, pks: Vec<ark_bls12_377::G1Affine>) -> Self;

    fn get_1st_round_register_polynomials(&self) -> Self::P1;
    //TODO: remove bitmask arg
    fn get_2nd_round_register_polynomials(&mut self, bitmask: Vec<Fr>, verifier_challenge: Fr) -> Self::P2;


    // TODO: move zeta_minus_omega_inv param to evaluations
    fn evaluate_register_polynomials(&self, point: Fr) -> E;
    // TODO: move zeta_minus_omega_inv param to evaluations
    fn compute_linearization_polynomial(&self, evaluations: &E, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr>;
    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>>;
    fn get_all_register_polynomials(self) -> Vec<DensePolynomial<Fr>>;

    //TODO: remove domains param
    fn compute_quotient_polynomial(&self, phi: Fr, domains: &Domains) -> DensePolynomial<Fr> {
        let w = utils::randomize(phi, &self.compute_constraint_polynomials());
        let (q_poly, r) = domains.compute_quotient(&w);
        assert_eq!(r, DensePolynomial::zero());
        q_poly
    }
}

pub trait RegisterPolynomials<E> {
    fn to_vec(self) -> Vec<DensePolynomial<Fr>>;
    fn evaluate(&self, point: Fr) -> E;
}

pub struct PackedAccountabilityRegisterPolynomials {
    pub c_poly: DensePolynomial<Fr>,
    pub acc_poly: DensePolynomial<Fr>,
}

impl PackedAccountabilityRegisterPolynomials {
    pub fn new(c_poly: DensePolynomial<Fr>, acc_poly: DensePolynomial<Fr>) -> Self {
        PackedAccountabilityRegisterPolynomials { c_poly, acc_poly }
    }
}

impl RegisterPolys for PackedAccountabilityRegisterPolynomials {
    type C = PackedRegisterCommitments;

    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> Self::C {
            PackedRegisterCommitments::new(
                f(&self.c_poly),
                f(&self.acc_poly),
            )
    }
}