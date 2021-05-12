use crate::piop::affine_addition::{AffineAdditionRegisters, AffineAdditionPolynomials, PartialSumsPolynomials, PartialSumsAndBitmaskPolynomials, AffineAdditionEvaluations, PartialSumsAndBitmaskCommitments};
use crate::piop::{ProverProtocol, RegisterPolynomials, RegisterEvaluations, RegisterCommitments, VerifierProtocol};
use crate::domains::Domains;
use ark_poly::polynomial::univariate::DensePolynomial;
use crate::{Bitmask, utils};
use ark_bw6_761::{Fr, G1Projective};
use crate::piop::bit_counting::BitCountingRegisters;

use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};
use crate::utils::{lagrange_evaluations, LagrangeEvaluations};
use ark_poly::Polynomial;
use ark_bls12_377::G1Affine;
use ark_ec::AffineCurve;


#[derive(CanonicalSerialize, CanonicalDeserialize)]
struct CountingCommitments {
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



struct CountingPolynomials {
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
struct CountingEvaluations {
    affine_addition_evaluations: AffineAdditionEvaluations,
    partial_counts_evaluation: Fr,
}

impl RegisterEvaluations for CountingEvaluations {
    fn as_vec(&self) -> Vec<Fr> {
        let mut evals = self.affine_addition_evaluations.as_vec();
        evals.push(self.partial_counts_evaluation);
        evals
    }
}



struct CountingScheme {
    affine_addition_registers: AffineAdditionRegisters,
    bit_counting_registers: BitCountingRegisters,
    register_evaluations: Option<CountingEvaluations>,
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
        CountingPolynomials {
            affine_addition_polynomials: self.affine_addition_registers.get_partial_sums_and_bitmask_polynomials(),
            partial_counts_polynomial: self.bit_counting_registers.get_partial_counts_polynomial()
        }
    }

    fn get_register_polynomials_to_commit2(&mut self, verifier_challenge: Fr) -> Self::P2 {
        ()
    }

    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<Fr>> {
        let affine_addition_polynomials = self.affine_addition_registers.get_register_polynomials();
        let partial_counts_polynomial = self.bit_counting_registers.get_partial_counts_polynomial();
        let mut polynomials = vec![];
        polynomials.extend(affine_addition_polynomials.to_vec());
        polynomials.push(partial_counts_polynomial);
        polynomials
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>> {
        let affine_addition_constraints = self.affine_addition_registers.compute_constraint_polynomials();
        let bit_counting_constraint = self.bit_counting_registers.compute_bit_counting_constraint();
        let mut constraints = vec![];
        constraints.extend(affine_addition_constraints);
        constraints.push(bit_counting_constraint);
        constraints
    }

    fn evaluate_register_polynomials(&mut self, point: Fr) -> Self::E {
        let affine_addition_evaluations = self.affine_addition_registers.evaluate_register_polynomials(point);
        let partial_counts_evaluation = self.bit_counting_registers.get_partial_counts_polynomial().evaluate(&point);
        let evals = CountingEvaluations {
            affine_addition_evaluations,
            partial_counts_evaluation,
        };
        self.register_evaluations = Some(evals.clone());
        evals
    }

    fn compute_linearization_polynomial(&self, phi: Fr, zeta: Fr) -> DensePolynomial<Fr> {
        let evals = self.register_evaluations.as_ref().unwrap();
        let affine_addition_parts =
            self.affine_addition_registers.compute_constraints_linearized(&evals.affine_addition_evaluations, zeta);
        let bit_counting_part =
            self.bit_counting_registers.compute_bit_counting_constraint_linearized();

        let mut parts = vec![];
        parts.extend(affine_addition_parts);
        parts.push(bit_counting_part);
        utils::randomize(phi, &parts)
    }
}



impl VerifierProtocol for CountingEvaluations {
    type C1 = CountingCommitments;
    type C2 = ();

    fn restore_commitment_to_linearization_polynomial(&self, phi: Fr, zeta_minus_omega_inv: Fr, commitments: &Self::C1, extra_commitments: &Self::C2) -> G1Projective {
        let powers_of_phi = utils::powers(phi, 6);
        let partial_sums_commitments= &commitments.affine_addition_commitments.partial_sums;
        let mut r_comm = self.affine_addition_evaluations.restore_commitment_to_linearization_polynomial(phi, zeta_minus_omega_inv, partial_sums_commitments, &());
        r_comm += commitments.partial_counts_commitment.mul(powers_of_phi[5]);
        r_comm
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use crate::tests::{random_bits, random_pks};
    use crate::KZG_BW6;

    #[test]
    fn test_polynomial_ordering() {
        let rng = &mut test_rng();
        let n = 16;
        let m = n - 1;

        let kzg_params = KZG_BW6::setup(m, rng);
        let mut scheme = CountingScheme::init(
            Domains::new(n),
            Bitmask::from_bits(&random_bits(m, 0.5, rng)),
            random_pks(m, rng),
        );

        let zeta = Fr::rand(rng);

        let actual_commitments = scheme.get_register_polynomials_to_commit1()
            .commit(|p| KZG_BW6::commit(&kzg_params.get_pk(), &p)).as_vec();
        let actual_evaluations = scheme.evaluate_register_polynomials(zeta).as_vec();
        let polynomials = scheme.get_register_polynomials_to_open();

        let expected_evaluations = polynomials.iter()
            .map(|p| p.evaluate(&zeta))
            .collect::<Vec<_>>();
        assert_eq!(actual_evaluations, expected_evaluations);


        let expected_commitments = polynomials.iter()
            .skip(2) // keyset commitment is publicly known
            .map(|p| KZG_BW6::commit(&kzg_params.get_pk(), &p))
            .collect::<Vec<_>>();
        assert_eq!(actual_commitments, expected_commitments);
    }
}