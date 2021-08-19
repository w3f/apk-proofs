use crate::piop::{ProverProtocol, RegisterEvaluations};
use ark_bw6_761::Fr;
use ark_poly::univariate::DensePolynomial;
use crate::domains::Domains;
use crate::{Bitmask, utils, AccountablePublicInput};
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
    type E = AffineAdditionEvaluationsWithoutBitmask;
    type PI = AccountablePublicInput;

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

    fn get_register_polynomials_to_commit2(&mut self, _verifier_challenge: Fr) -> () {
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

    // bitmask register polynomial is not committed to...
    fn evaluate_register_polynomials(&mut self, point: Fr) -> AffineAdditionEvaluationsWithoutBitmask {
        let evals: AffineAdditionEvaluations = self.registers.evaluate_register_polynomials(point);
        self.register_evaluations = Some(evals.clone());
        AffineAdditionEvaluationsWithoutBitmask {
            keyset: evals.keyset,
            partial_sums: evals.partial_sums,
        }
    }

    fn compute_linearization_polynomial(&self, phi: Fr, zeta: Fr) -> DensePolynomial<Fr> {
        let evals = self.register_evaluations.as_ref().unwrap();
        let parts = self.registers.compute_constraints_linearized(evals, zeta);
        utils::randomize(phi, &parts)
    }
}