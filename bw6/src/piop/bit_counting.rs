use ark_poly::{Evaluations, Radix2EvaluationDomain, EvaluationDomain};
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
// let L[i] = 0, i = 0,...,n-2, and L[n-1] = S.
// The constraint then is
// r[i+1] = r[i] + b[i] - L[i], i = 0,...,n-1
// For i = 0,...,n-2 it ensures that the register is well-formed: r[i+1] = r[i] + b[i],
// for i = n-1, r[n] = r[n-1] + b[n-1] - S, but r[n] = r[0] and r[n-1] + b[n-1] = S.
// It follows that r[0] = 0 and S is properly constrained.

// After interpolation L this defined becomes S * L_{n-1}(Z) that is known by the verifier, and easy to compute.

pub(crate) struct BitCountingRegisters {
    domain: Radix2EvaluationDomain<Fr>,

    bitmask: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    partial_counts: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    partial_counts_shifted: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
}

impl BitCountingRegisters {
    pub fn new(bitmask: &Bitmask) -> Self {
        let bitmask = bitmask.to_bits_as_field_elements();
        let partial_counts = Self::build_partial_counts_register(&bitmask);
        let mut partial_counts_shifted = partial_counts.clone();
        partial_counts_shifted.rotate_left(1);
        Self::new_unchecked(bitmask, partial_counts, partial_counts_shifted)
    }

    fn new_unchecked(bitmask: Vec<Fr>,
                     partial_counts: Vec<Fr>,
                     partial_counts_shifted: Vec<Fr>,
    ) -> Self {
        let domain = Radix2EvaluationDomain::<Fr>::new(bitmask.len()).unwrap();
        Self {
            domain,

            bitmask: Evaluations::from_vec_and_domain(bitmask, domain),
            partial_counts: Evaluations::from_vec_and_domain(partial_counts, domain),
            partial_counts_shifted: Evaluations::from_vec_and_domain(partial_counts_shifted, domain),
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
        let count = self.bitmask.evals.iter().sum();
        let n = self.domain.size();
        let mut count_at_l_last = vec![Fr::zero(); n];
        count_at_l_last[n-1] = count;
        let count_at_l_last = Evaluations::from_vec_and_domain(count_at_l_last, self.domain);
        let constraint = &(&(&self.partial_counts_shifted - &self.partial_counts) - &self.bitmask) + &count_at_l_last;
        constraint.interpolate()
    }

    pub fn compute_bit_counting_constraint_linearized(&self) -> DensePolynomial<Fr> {
        self.get_partial_counts_polynomial()
    }

    pub fn get_partial_counts_polynomial(&self) -> DensePolynomial<Fr> {
        self.partial_counts.interpolate_by_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;
    use crate::tests::{random_bitmask, random_bits};
    use ark_poly::Polynomial;

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
        // assert_eq!(constraint.degree(), n - 1);
        assert!(constraint.divide_by_vanishing_poly(registers.domain).unwrap().1.is_zero());
    }
}