use crate::piop::ProverProtocol;
use ark_bw6_761::Fr;
use ark_poly::univariate::DensePolynomial;
use crate::domains::Domains;
use crate::Bitmask;
use crate::piop::affine_addition::{AffineAdditionRegisters, AffineAdditionEvaluations, PartialSumsPolynomials};

pub struct BasicRegisterBuilder {
    registers: AffineAdditionRegisters,
    register_evaluations: Option<AffineAdditionEvaluations>,
}

impl ProverProtocol for BasicRegisterBuilder {
    type P1 = PartialSumsPolynomials;
    type P2 = ();
    type E = AffineAdditionEvaluations;

    fn init(domains: Domains, bitmask: Bitmask, pks: Vec<ark_bls12_377::G1Affine>) -> Self {
        BasicRegisterBuilder {
            registers:  AffineAdditionRegisters::new(domains, &bitmask, pks),
            register_evaluations: None,
        }
    }

    fn get_register_polynomials_to_commit(&self) -> PartialSumsPolynomials {
        self.registers.get_partial_sums_register_polynomials()
    }

    fn get_register_polynomials_to_commit_extra(&mut self, verifier_challenge: Fr) -> () {
        ()
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>> {
        self.registers.compute_constraint_polynomials()
    }

    fn evaluate_register_polynomials(&mut self, point: Fr) -> AffineAdditionEvaluations {
        let mut evals = self.registers.evaluate_register_polynomials(point);
        self.register_evaluations = Some(evals.clone());
        evals.bitmask = None;
        evals
    }

    fn compute_linearization_polynomial(&self, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr> {
        self.registers.compute_linearization_polynomial(self.register_evaluations.as_ref().unwrap(), phi, zeta_minus_omega_inv)
    }

    fn get_all_register_polynomials(self) -> Vec<DensePolynomial<Fr>> {
        self.registers.get_register_polynomials().to_vec()
    }
}