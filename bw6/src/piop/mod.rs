use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_bw6_761::Fr;

use crate::domains::Domains;
use crate::utils;
use crate::constraints::Registers;

pub trait Piop<E> {
    // TODO: move zeta_minus_omega_inv param to evaluations
    fn evaluate_register_polynomials(&self, point: Fr) -> E;
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

pub trait PiopDecorator<E>: Piop<E> {
    // TODO: move zeta_minus_omega_inv param to evaluations
    fn wrap(registers: Registers, bitmask: Vec<Fr>, bitmask_chunks_aggregation_challenge: Fr) -> Self;
    fn get_accountable_register_polynomials(&self) -> Option<PackedAccountabilityRegisterPolynomials>;
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