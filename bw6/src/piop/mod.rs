use ark_bw6_761::{Fr, G1Affine};
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use crate::{Bitmask, utils};
use crate::domains::Domains;
use crate::utils::LagrangeEvaluations;

pub mod packed;
pub mod affine_addition;
pub mod basic;
pub mod bitmask_packing;

pub trait RegisterCommitments: CanonicalSerialize + CanonicalDeserialize {
    fn as_vec(&self) -> Vec<G1Affine>;
}

pub trait RegisterPolynomials {
    type C: RegisterCommitments;
    fn commit<F: Clone + Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> Self::C;
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

// Represents a polynomial protocol as seen by the prover.
pub trait ProverProtocol {
    type P1: RegisterPolynomials;
    type P2: RegisterPolynomials;
    type E: RegisterEvaluations;

    fn init(domains: Domains, bitmask: Bitmask, pks: Vec<ark_bls12_377::G1Affine>) -> Self;

    // These 2 methods together return register polynomials the prover should commit to.
    // The 2nd one is used only in the "packed" scheme as it requires an additional challenge
    // (to aggregate the bitmask chunks) from the verifier,
    // that can be received only after the bitmask has been committed.
    fn get_register_polynomials_to_commit1(&self) -> Self::P1;
    fn get_register_polynomials_to_commit2(&mut self, verifier_challenge: Fr) -> Self::P2;

    // This method returns register polynomials the prover should open. Those are the same polynomials
    // as the previous 2 methods together, and additionally 2 polynomials representing the keyset
    // (prover doesn't need to commit to them, as verifier knows them anyway, but still should open).
    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<Fr>>;

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>>;

    //TODO: remove domains param
    fn compute_quotient_polynomial(&self, phi: Fr, domains: &Domains) -> DensePolynomial<Fr> {
        let w = utils::randomize(phi, &self.compute_constraint_polynomials());
        let (q_poly, r) = domains.compute_quotient(&w);
        assert_eq!(r, DensePolynomial::zero());
        q_poly
    }

    fn evaluate_register_polynomials(&mut self, point: Fr) -> Self::E;

    // TODO: move zeta_minus_omega_inv param to evaluations
    fn compute_linearization_polynomial(&self, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr>;
}

pub trait RegisterEvaluations: CanonicalSerialize + CanonicalDeserialize {
    fn as_vec(&self) -> Vec<Fr>;
}

pub trait VerifierProtocol {
    type AC: RegisterCommitments;
    type C: RegisterCommitments;

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
}
