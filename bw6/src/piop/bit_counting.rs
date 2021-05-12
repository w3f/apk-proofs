use ark_poly::{Evaluations, Radix2EvaluationDomain, EvaluationDomain, Polynomial};
use ark_bw6_761::Fr;
use crate::Bitmask;
use ark_ff::Zero;
use ark_std::iter::once;
use ark_poly::univariate::DensePolynomial;


// This "gadget" is used in the 'counting' scheme to constraint the number of set bits in the bitmask.

// It adds a single register tracking partial bitmask sums:
// r[0] = 0, r[i] = b[0] + ... + b[i-1], i = 1,...,n-1,
// and a single constraint operating on this register and the bitmask register b.
// Let S = b[0]+...+b[n-1] - the value we want to constraint and
// let SL[i] = 0, i = 0,...,n-2, and SL[n-1] = S.
// The constraint then is
// r[i+1] = r[i] + b[i] - SL[i], i = 0,...,n-1
// For i = 0,...,n-2 it ensures that the register is well-formed: r[i+1] = r[i] + b[i],
// for i = n-1, r[n] = r[n-1] + b[n-1] - S, but r[n] = r[0] and r[n-1] + b[n-1] = S.
// It follows that r[0] = 0 and S is properly constrained.

// After interpolation SL this defined becomes S * L_{n-1}(Z) that is known by the verifier, and easy to compute.

// What makes this constraint different from any other in the codebase is that it is linear in the registers involved,
// so by construction deg(C) <= n-1, and C(w^i) = 0 for i = 0,...,n-1, that results in C = 0 identically.
// To check this constraint holds, the verifier doesn't need to check C(Z) = Q(Z)(Z^n - 1), it could just check C(zeta) = 0,
// But we use the former check not to handle this case differently.

pub(crate) struct BitCountingRegisters {
    domain: Radix2EvaluationDomain<Fr>,

    bitmask: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    partial_counts: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
}

impl BitCountingRegisters {
    pub fn new(bitmask: &Bitmask) -> Self {
        let bitmask = bitmask.to_bits_as_field_elements();
        let partial_counts = Self::build_partial_counts_register(&bitmask);
        Self::new_unchecked(bitmask, partial_counts)
    }

    fn new_unchecked(bitmask: Vec<Fr>,
                     partial_counts: Vec<Fr>,
    ) -> Self {
        let domain = Radix2EvaluationDomain::<Fr>::new(bitmask.len()).unwrap();
        Self {
            domain,
            bitmask: Evaluations::from_vec_and_domain(bitmask, domain),
            partial_counts: Evaluations::from_vec_and_domain(partial_counts, domain),
        }
    }

    /// Returns length n vec (0, b[0], b[0] + b[1], ..., b[0] + b[1] + b[n-2])
    fn build_partial_counts_register(bitmask: &[Fr]) -> Vec<Fr> {
        let partial_counts = bitmask.iter().scan(Fr::zero(), |state, bit| {
            *state += bit;
            Some(*state)
        });
        once(Fr::zero())
            .chain(partial_counts.take(bitmask.len() - 1))
            .collect()
    }

    pub fn compute_bit_counting_constraint(&self) -> DensePolynomial<Fr> {
        DensePolynomial::zero()
    }

    // Though the constraint is zero, the verifier still needs the opening of the register in zeta * omega.
    pub fn compute_bit_counting_constraint_linearized(&self) -> DensePolynomial<Fr> {
        self.get_partial_counts_polynomial()
    }

    pub fn get_partial_counts_polynomial(&self) -> DensePolynomial<Fr> {
        //TODO: cache
        self.partial_counts.interpolate_by_ref()
    }

    pub fn evaluate_register(&self, zeta: Fr) -> BitCountingEvaluation {
        let eval = self.get_partial_counts_polynomial().evaluate(&zeta);
        BitCountingEvaluation(eval)
    }

    #[cfg(test)]
    fn get_bitmask_polynomial(&self) -> DensePolynomial<Fr> {
        self.bitmask.interpolate_by_ref()
    }
}

pub struct BitCountingEvaluation(Fr);

impl BitCountingEvaluation {
    fn evaluate_constraint(
        &self,
        count: Fr,
        partial_counts_at_zeta_omega: Fr,
        bitmask_at_zeta: Fr,
        l_last_at_zeta: Fr,
    ) -> Fr {
        partial_counts_at_zeta_omega - self.0 - bitmask_at_zeta + count * l_last_at_zeta
    }

    fn evaluate_constraint_at_zeta(
        &self,
        count: Fr,
        bitmask_at_zeta: Fr,
        l_last_at_zeta: Fr,
    ) -> Fr {
        self.evaluate_constraint(count, Fr::zero(), bitmask_at_zeta, l_last_at_zeta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use crate::tests::{random_bitmask, random_bits};
    use ark_poly::Polynomial;
    use crate::utils;

    #[test]
    fn test_partial_counts_register() {
        let rng = &mut test_rng();
        let n = 16;
        let bitmask = random_bitmask(n, rng);
        let partial_counts = BitCountingRegisters::build_partial_counts_register(&bitmask);

        assert_eq!(partial_counts.len(), 16);
        assert_eq!(partial_counts[0], Fr::zero());
        for i in 1..n-1 {
            assert_eq!(partial_counts[i], bitmask[..i].iter().sum());
        }
    }

    #[test]
    fn test_bit_counting_constraint() {
        let rng = &mut test_rng();
        let n = 16;
        let bitmask = Bitmask::from_bits(&random_bits(n, 2.0 / 3.0, rng));
        let registers = BitCountingRegisters::new(&bitmask);

        let constraint = registers.compute_bit_counting_constraint();
        assert!(constraint.divide_by_vanishing_poly(registers.domain).unwrap().1.is_zero());
    }

    #[test]
    fn test_bit_counting_constraint_eval() {
        let rng = &mut test_rng();
        let n = 16;
        let domain = Radix2EvaluationDomain::<Fr>::new(n).unwrap();
        let bitmask = Bitmask::from_bits(&random_bits(n, 2.0 / 3.0, rng));
        let registers = BitCountingRegisters::new(&bitmask);

        let zeta = Fr::rand(rng);

        let constraint = registers.compute_bit_counting_constraint();
        let expected_constraint_eval = constraint.evaluate(&zeta);

        let omega = domain.group_gen;
        let count = Fr::from(bitmask.count_ones() as u8);
        let evals_at_zeta = utils::lagrange_evaluations(zeta, domain);
        let partial_counts_at_zeta_omega = registers.get_partial_counts_polynomial().evaluate(&(zeta * omega));
        let bitmask_at_zeta = registers.get_bitmask_polynomial().evaluate(&zeta);

        let register_eval = registers.evaluate_register(zeta);
        let actual_constraint_eval = register_eval.evaluate_constraint(count, partial_counts_at_zeta_omega, bitmask_at_zeta, evals_at_zeta.l_last);
        assert_eq!(actual_constraint_eval, expected_constraint_eval);

        let constraint_at_zeta = register_eval.evaluate_constraint_at_zeta(count, bitmask_at_zeta, evals_at_zeta.l_last);
        let constraint_at_zeta_omega = registers.compute_bit_counting_constraint_linearized().evaluate(&(zeta * omega));
        assert_eq!(constraint_at_zeta + constraint_at_zeta_omega, expected_constraint_eval);
    }
}