use crate::piop::affine_addition::{AffineAdditionRegisters, AffineAdditionPolynomials};
use crate::piop::{ProverProtocol, RegisterPolynomials, RegisterEvaluations, RegisterCommitments};
use crate::domains::Domains;
use ark_poly::polynomial::univariate::DensePolynomial;
use crate::Bitmask;
use ark_bw6_761::{Fr};
use crate::piop::bit_counting::{BitCountingRegisters, BitCountingPolynomials};

use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};


struct CountingScheme {
    affine_addition_registers: AffineAdditionRegisters,
    bit_counting_registers: BitCountingRegisters,
    register_evaluations: Option<CountingEvaluations>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
struct CountingCommitments {

}

impl RegisterCommitments for CountingCommitments {
    fn as_vec(&self) -> Vec<ark_bw6_761::G1Affine> {
        unimplemented!()
    }
}

struct CountingPolynomials {
    affine_addition_polynomials: AffineAdditionPolynomials,
    bit_counting_polynomials: BitCountingPolynomials,
}

impl RegisterPolynomials for CountingPolynomials {
    type C = CountingCommitments;

    fn commit<F: Clone + Fn(&DensePolynomial<Fr>) -> ark_bw6_761::G1Affine>(&self, f: F) -> Self::C {
        unimplemented!()
    }
}

impl ProverProtocol for CountingScheme {
    type P1 = CountingPolynomials;
    type P2 = ();
    type E = CountingEvaluations;

    fn init(domains: Domains, bitmask: Bitmask, pks: Vec<ark_bls12_377::G1Affine>) -> Self {
        CountingScheme {
            affine_addition_registers: AffineAdditionRegisters::new(domains, &bitmask, pks),
            bit_counting_registers: BitCountingRegisters::new(&bitmask),
            register_evaluations: None,
        }
    }

    fn get_register_polynomials_to_commit1(&self) -> Self::P1 {
        unimplemented!()
    }

    fn get_register_polynomials_to_commit2(&mut self, verifier_challenge: Fr) -> Self::P2 {
        unimplemented!()
    }

    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<Fr>> {
        unimplemented!()
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>> {
        unimplemented!()
    }

    fn evaluate_register_polynomials(&mut self, point: Fr) -> Self::E {
        unimplemented!()
    }

    fn compute_linearization_polynomial(&self, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr> {
        unimplemented!()
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
struct CountingEvaluations {

}

impl RegisterEvaluations for CountingEvaluations {
    fn as_vec(&self) -> Vec<Fr> {
        unimplemented!()
    }
}