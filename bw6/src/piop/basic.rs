use crate::piop::Protocol;
use ark_bw6_761::Fr;
use ark_poly::univariate::DensePolynomial;
use crate::domains::Domains;
use crate::Bitmask;
use crate::piop::affine_addition::{AffineAdditionRegisters, AffineAdditionEvaluations, PartialSumsPolynomials};

pub struct BasicRegisterBuilder {
    affine_addition_registers: AffineAdditionRegisters,
}

impl Protocol for BasicRegisterBuilder {
    type P1 = PartialSumsPolynomials;
    type P2 = ();
    type E = AffineAdditionEvaluations;

    fn init(domains: Domains, bitmask: &Bitmask, pks: Vec<ark_bls12_377::G1Affine>) -> Self {
        BasicRegisterBuilder {
            affine_addition_registers:  AffineAdditionRegisters::new(domains, bitmask, pks)
        }
    }

    fn get_1st_round_register_polynomials(&self) -> PartialSumsPolynomials {
        self.affine_addition_registers.get_partial_sums_register_polynomials()
    }

    fn get_2nd_round_register_polynomials(&mut self, bitmask: Vec<Fr>, verifier_challenge: Fr) -> () {
        ()
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>> {
        self.affine_addition_registers.compute_constraint_polynomials()
    }

    fn evaluate_register_polynomials(&self, point: Fr) -> AffineAdditionEvaluations {
        self.affine_addition_registers.evaluate_register_polynomials(point)
    }

    fn compute_linearization_polynomial(&self, evaluations: &AffineAdditionEvaluations, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr> {
        self.affine_addition_registers.compute_linearization_polynomial(evaluations, phi, zeta_minus_omega_inv)
    }

    fn get_all_register_polynomials(self) -> Vec<DensePolynomial<Fr>> {
        self.affine_addition_registers.get_all_register_polynomials()
    }
}