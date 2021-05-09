use crate::piop::{ProverProtocol, RegisterEvaluations};
use ark_bw6_761::Fr;
use ark_poly::univariate::DensePolynomial;
use crate::domains::Domains;
use crate::Bitmask;
use crate::piop::affine_addition::{AffineAdditionRegisters, AffineAdditionEvaluations, PartialSumsPolynomials};

use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct AffineAdditionEvaluationsWithoutBitmask {
    pub keyset: (Fr, Fr),
    pub partial_sums: (Fr, Fr),
}

impl RegisterEvaluations for AffineAdditionEvaluationsWithoutBitmask {
    fn as_vec(&self) -> Vec<Fr> {
        vec![
            self.keyset.0,
            self.keyset.1,
            self.partial_sums.0,
            self.partial_sums.1,
        ]
    }
}

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

    fn get_register_polynomials_to_commit1(&self) -> PartialSumsPolynomials {
        let polys = self.registers.get_register_polynomials();
        PartialSumsPolynomials(
            polys.partial_sums.0,
            polys.partial_sums.1,
        )
    }

    fn get_register_polynomials_to_commit2(&mut self, verifier_challenge: Fr) -> () {
        ()
    }

    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<Fr>> {
        let polys = self.registers.get_register_polynomials();
        vec![
            polys.keyset.0,
            polys.keyset.1,
            polys.partial_sums.0,
            polys.partial_sums.1,
        ]
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>> {
        self.registers.compute_constraint_polynomials()
    }

    fn evaluate_register_polynomials(&mut self, point: Fr) -> AffineAdditionEvaluations {
        let mut evals = self.registers.evaluate_register_polynomials(point);
        self.register_evaluations = Some(evals.clone());
        evals
    }

    fn compute_linearization_polynomial(&self, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr> {
        self.registers.compute_linearization_polynomial(self.register_evaluations.as_ref().unwrap(), phi, zeta_minus_omega_inv)
    }
}