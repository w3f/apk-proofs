use ark_poly::{Evaluations, Polynomial};
use ark_bw6_761::Fr;
use crate::Bitmask;
use ark_ff::{Zero, One};
use ark_std::iter::once;
use ark_poly::univariate::DensePolynomial;

use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};
use crate::domains::Domains;


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
    domains: Domains,
    bitmask: Vec<Fr>,
    partial_counts: DensePolynomial<Fr>,
}

impl BitCountingRegisters {
    pub fn new(domains: Domains, bitmask: &Bitmask) -> Self {
        let mut bitmask = bitmask.to_bits_as_field_elements();
        bitmask.resize_with(domains.size, || Fr::zero());
        let partial_counts = Self::build_partial_counts_register(&bitmask);
        Self::new_unchecked(domains, bitmask, partial_counts)
    }

    fn new_unchecked(domains: Domains,
                     bitmask: Vec<Fr>,
                     partial_counts: Vec<Fr>,
    ) -> Self {
        let partial_counts= Evaluations::from_vec_and_domain(partial_counts, domains.domain.clone()).interpolate();
        Self {
            domains,
            bitmask,
            partial_counts,
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

    pub fn constraints(&self) -> Vec<DensePolynomial<Fr>> {
        vec![
            BitCount::constraint_poly(),
            BitmaskEndsWithZero::constraint_poly(self),
        ]
    }

    pub fn constraints_lin(&self) -> Vec<DensePolynomial<Fr>> {
        vec![
            BitCount::linearization(&self),
            BitmaskEndsWithZero::linearization(),
        ]
    }

    pub fn get_partial_counts_polynomial(&self) -> DensePolynomial<Fr> {
        self.partial_counts.clone()
    }

    pub fn evaluate_partial_counts_register(&self, zeta: Fr) -> BitCountingEvaluation {
        let eval = self.partial_counts.evaluate(&zeta);
        BitCountingEvaluation(eval)
    }

    #[cfg(test)]
    fn get_bitmask_polynomial(&self) -> DensePolynomial<Fr> {
        Evaluations::from_vec_and_domain( self.bitmask.clone(), self.domains.domain.clone()).interpolate()
    }
}

struct BitmaskEndsWithZero {}

// Constraints the last bitmask element to 0
// by showing that the polynomial b(X) * L_{n-1}(X) is zero over the domain.
impl BitmaskEndsWithZero {

    // C = b * L_{n-1}
    fn constraint_poly(registers: &BitCountingRegisters) -> DensePolynomial<Fr> {
        let n = registers.domains.size;
        let mut ln = vec![Fr::zero(); n];
        ln[n-1] = Fr::one();
        let ln_x2 = registers.domains.amplify_x2(ln);
        let b_x2 = registers.domains.amplify_x2(registers.bitmask.clone());
        let c = &b_x2 * &ln_x2;
        c.interpolate()
    }

    fn linearization() -> DensePolynomial<Fr> {
        DensePolynomial::zero()
    }

    // C(z) = b(z) * L_{n-1}(z)
    fn _evaluate_full(bitmask_at_zeta: Fr, l_last_at_zeta: Fr) -> Fr {
        bitmask_at_zeta * l_last_at_zeta
    }

    fn evaluate_main(bitmask_at_zeta: Fr, l_last_at_zeta: Fr) -> Fr {
        Self::_evaluate_full(bitmask_at_zeta, l_last_at_zeta)
    }
}

struct BitCount {}

// partial_counts(wZ) - partial_counts(Z) - bitmask(Z) + count*L_{n-1}(Z)
impl BitCount {

    fn constraint_poly() -> DensePolynomial<Fr> {
        DensePolynomial::zero()
    }

    // Though the constraint is zero, the verifier still needs the opening of the register in zeta * omega.
    fn linearization(registers: &BitCountingRegisters) -> DensePolynomial<Fr> {
        registers.get_partial_counts_polynomial()
    }

    fn _evaluate_full(
        eval: &BitCountingEvaluation,
        count: Fr,
        partial_counts_at_zeta_omega: Fr,
        bitmask_at_zeta: Fr,
        l_last_at_zeta: Fr,
    ) -> Fr {
        partial_counts_at_zeta_omega - eval.0 - bitmask_at_zeta + count * l_last_at_zeta
    }

    fn evaluate_main(
        eval: &BitCountingEvaluation,
        count: Fr,
        bitmask_at_zeta: Fr,
        l_last_at_zeta: Fr,
    ) -> Fr {
        Self::_evaluate_full(eval, count, Fr::zero(), bitmask_at_zeta, l_last_at_zeta)
    }
}


#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct BitCountingEvaluation(pub Fr);

impl BitCountingEvaluation {
    pub fn evaluate_constraints_at_zeta(
        &self,
        count: Fr,
        bitmask_at_zeta: Fr,
        l_last_at_zeta: Fr,
    ) -> Vec<Fr> {
        vec![
            BitCount::evaluate_main(&self, count, bitmask_at_zeta, l_last_at_zeta),
            BitmaskEndsWithZero::evaluate_main(bitmask_at_zeta, l_last_at_zeta),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use crate::tests::{random_bitmask, random_bits};
    use ark_poly::Polynomial;
    use crate::utils;
    use crate::utils::lagrange_evaluations;

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
        let domains = Domains::new(n);

        let bitmask = Bitmask::from_bits(&random_bits(n, 2.0 / 3.0, rng));
        let count = Fr::from(bitmask.count_ones() as u8);
        let registers = BitCountingRegisters::new(domains.clone(), &bitmask);

        let z = Fr::rand(rng);
        let w = domains.omega;

        let acc_z = registers.evaluate_partial_counts_register(z);
        let acc_zw = registers.get_partial_counts_polynomial().evaluate(&(z * w));
        let bitmask_z = registers.get_bitmask_polynomial().evaluate(&z);
        let domain_z = utils::lagrange_evaluations(z, domains.domain);

        let x_full = BitCount::constraint_poly();
        let x_lin = BitCount::linearization(&registers);
        let x_eval_full = BitCount::_evaluate_full(&acc_z, count, acc_zw, bitmask_z, domain_z.l_last);
        let x_eval_main = BitCount::evaluate_main(&acc_z, count, bitmask_z, domain_z.l_last);
        let x_lin_zw = x_lin.evaluate(&(z * w));

        assert_eq!(x_eval_full, x_full.evaluate(&z));
        assert_eq!(x_eval_full, x_eval_main + x_lin_zw);
        assert!(x_full.divide_by_vanishing_poly(domains.domain).unwrap().1.is_zero()); // actually x_full is 0 over the field
    }

    #[test]
    fn test_bitmask_ends_with_zero_constraint() {
        let rng = &mut test_rng();
        let n = 16;
        let domains = Domains::new(n);
        let domain = domains.domain;

        let bits = random_bits(n, 2.0 / 3.0, rng);

        let mut good_bitmask = bits.clone();
        good_bitmask[n-1] = false;
        let good_bitmask = Bitmask::from_bits(&good_bitmask);
        let registers = BitCountingRegisters::new(domains.clone(), &good_bitmask);
        let constraint = BitmaskEndsWithZero::constraint_poly(&registers);
        assert_eq!(constraint.degree(), 2 * (n - 1));
        assert!(constraint.divide_by_vanishing_poly(domain).unwrap().1.is_zero());

        let zeta = Fr::rand(rng);
        let prover_eval = constraint.evaluate(&zeta);
        let domain_evals = lagrange_evaluations(zeta, domain);
        let bitmask_at_zeta = registers.get_bitmask_polynomial().evaluate(&zeta);
        let verifier_eval = BitmaskEndsWithZero::_evaluate_full(bitmask_at_zeta, domain_evals.l_last);
        assert_eq!(prover_eval, verifier_eval);

        let mut bad_bitmask = bits.clone();
        bad_bitmask[n-1] = true;
        let bad_bitmask = Bitmask::from_bits(&bad_bitmask);
        let registers = BitCountingRegisters::new(domains.clone(), &bad_bitmask);
        let constraint = BitmaskEndsWithZero::constraint_poly(&registers);
        assert!(!constraint.divide_by_vanishing_poly(domain).unwrap().1.is_zero());
    }
}