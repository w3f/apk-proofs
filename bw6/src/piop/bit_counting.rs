use ark_poly::{Evaluations, Radix2EvaluationDomain, EvaluationDomain};
use ark_bw6_761::Fr;
use crate::Bitmask;
use ark_ff::Zero;
use ark_std::iter::once;


pub(crate) struct BitCountingRegisters {
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
}

pub(crate) struct BitCountingPolynomials {}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use crate::tests::random_bitmask;

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
}