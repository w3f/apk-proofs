use ark_bw6_761::{Fr, G1Projective};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use crate::{Bitmask, CountingPublicInput, Keyset, utils};
use crate::domains::Domains;
use crate::piop::{ProverProtocol, RegisterCommitments, RegisterEvaluations, RegisterPolynomials, VerifierProtocol};
use crate::piop::affine_addition::{AffineAdditionEvaluations, AffineAdditionRegisters, PartialSumsAndBitmaskCommitments, PartialSumsAndBitmaskPolynomials};
use crate::piop::bit_counting::{BitCountingEvaluation, BitCountingRegisters};
use crate::utils::LagrangeEvaluations;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct CountingCommitments {
    affine_addition_commitments: PartialSumsAndBitmaskCommitments,
    partial_counts_commitment: ark_bw6_761::G1Affine,
}

impl RegisterCommitments for CountingCommitments {
    fn as_vec(&self) -> Vec<ark_bw6_761::G1Affine> {
        let mut commitments = self.affine_addition_commitments.as_vec();
        commitments.push(self.partial_counts_commitment);
        commitments
    }
}


pub struct CountingPolynomials {
    affine_addition_polynomials: PartialSumsAndBitmaskPolynomials,
    partial_counts_polynomial: DensePolynomial<Fr>,
}

impl RegisterPolynomials for CountingPolynomials {
    type C = CountingCommitments;

    fn commit<F: Clone + Fn(&DensePolynomial<Fr>) -> ark_bw6_761::G1Affine>(&self, f: F) -> Self::C {
        CountingCommitments {
            affine_addition_commitments: self.affine_addition_polynomials.commit(f.clone()),
            partial_counts_commitment: f(&self.partial_counts_polynomial),
        }
    }
}


#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct CountingEvaluations {
    affine_addition_evaluations: AffineAdditionEvaluations,
    partial_counts_evaluation: BitCountingEvaluation,
}

impl RegisterEvaluations for CountingEvaluations {
    fn as_vec(&self) -> Vec<Fr> {
        let mut evals = self.affine_addition_evaluations.as_vec();
        evals.push(self.partial_counts_evaluation.0);
        evals
    }
}


pub struct CountingScheme {
    affine_addition_registers: AffineAdditionRegisters,
    bit_counting_registers: BitCountingRegisters,
    register_evaluations: Option<CountingEvaluations>,
}

impl ProverProtocol for CountingScheme {
    type P1 = CountingPolynomials;
    type P2 = ();
    type E = CountingEvaluations;
    type PI = CountingPublicInput;

    fn init(domains: Domains, bitmask: Bitmask, keyset: Keyset) -> Self {
        CountingScheme {
            affine_addition_registers: AffineAdditionRegisters::new(domains.clone(), keyset, &bitmask.to_bits()),
            bit_counting_registers: BitCountingRegisters::new(domains, &bitmask),
            register_evaluations: None,
        }
    }

    fn get_register_polynomials_to_commit1(&self) -> Self::P1 {
        CountingPolynomials {
            affine_addition_polynomials: self.affine_addition_registers.get_partial_sums_and_bitmask_polynomials(),
            partial_counts_polynomial: self.bit_counting_registers.get_partial_counts_polynomial(),
        }
    }

    fn get_register_polynomials_to_commit2(&mut self, _verifier_challenge: Fr) -> Self::P2 {
        ()
    }

    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<Fr>> {
        [
            self.affine_addition_registers.get_register_polynomials().to_vec(),
            vec![self.bit_counting_registers.get_partial_counts_polynomial()],
        ].concat()
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>> {
        [
            self.affine_addition_registers.compute_constraint_polynomials(),
            self.bit_counting_registers.constraints(),
        ].concat()
    }

    fn evaluate_register_polynomials(&mut self, point: Fr) -> Self::E {
        let affine_addition_evaluations = self.affine_addition_registers.evaluate_register_polynomials(point);
        let partial_counts_evaluation = self.bit_counting_registers.evaluate_partial_counts_register(point);
        let evals = CountingEvaluations {
            affine_addition_evaluations,
            partial_counts_evaluation,
        };
        self.register_evaluations = Some(evals.clone());
        evals
    }

    fn compute_linearization_polynomial(&self, phi: Fr, zeta: Fr) -> DensePolynomial<Fr> {
        let evals = self.register_evaluations.as_ref().unwrap();
        let parts = [
            self.affine_addition_registers.compute_constraints_linearized(&evals.affine_addition_evaluations, zeta),
            self.bit_counting_registers.constraints_lin(),
        ].concat();
        utils::randomize(phi, &parts)
    }
}


impl VerifierProtocol for CountingEvaluations {
    type C1 = CountingCommitments;
    type C2 = ();

    const POLYS_OPENED_AT_ZETA: usize = 7;

    fn restore_commitment_to_linearization_polynomial(&self, phi: Fr, zeta_minus_omega_inv: Fr, commitments: &Self::C1, _extra_commitments: &Self::C2) -> G1Projective {
        let powers_of_phi = utils::powers(phi, 6);
        let partial_sums_commitments = &commitments.affine_addition_commitments.partial_sums;
        let mut r_comm = self.affine_addition_evaluations.restore_commitment_to_linearization_polynomial(phi, zeta_minus_omega_inv, partial_sums_commitments, &());
        r_comm += commitments.partial_counts_commitment * powers_of_phi[5];
        r_comm
    }
}


impl CountingEvaluations {
    pub fn evaluate_constraint_polynomials(&self,
                                           apk: ark_bls12_377::G1Affine,
                                           count: Fr,
                                           evals_at_zeta: &LagrangeEvaluations<Fr>,
    ) -> Vec<Fr> {
        let b_at_zeta = self.affine_addition_evaluations.bitmask;
        [
            self.affine_addition_evaluations.evaluate_constraint_polynomials(apk, evals_at_zeta),
            self.partial_counts_evaluation.evaluate_constraints_at_zeta(count, b_at_zeta, evals_at_zeta.l_last),
        ].concat()
    }
}


#[cfg(test)]
mod tests {
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};
    use fflonk::pcs::{PCS, PcsParams};

    use crate::NewKzgBw6;
    use crate::test_helpers::{_random_bits, random_pks};

    use super::*;

    #[test]
    fn test_polynomial_ordering() {
        let rng = &mut test_rng();
        let n = 16;
        let m = n - 1;


        let kzg_params = NewKzgBw6::setup(m, rng);
        let mut keyset = Keyset::new(random_pks(m, rng));
        keyset.amplify();
        let mut scheme = CountingScheme::init(
            Domains::new(n),
            Bitmask::from_bits(&_random_bits(m, 0.5, rng)),
            keyset,
        );

        let zeta = Fr::rand(rng);

        let actual_commitments = scheme.get_register_polynomials_to_commit1()
            .commit(|p| NewKzgBw6::commit(&kzg_params.ck(), &p).0).as_vec();
        let actual_evaluations = scheme.evaluate_register_polynomials(zeta).as_vec();
        let polynomials = scheme.get_register_polynomials_to_open();

        let expected_evaluations = polynomials.iter()
            .map(|p| p.evaluate(&zeta))
            .collect::<Vec<_>>();
        assert_eq!(actual_evaluations, expected_evaluations);


        let expected_commitments = polynomials.iter()
            .skip(2) // keyset commitment is publicly known
            .map(|p| NewKzgBw6::commit(&kzg_params.ck(), &p).0)
            .collect::<Vec<_>>();
        assert_eq!(actual_commitments, expected_commitments);
    }
}