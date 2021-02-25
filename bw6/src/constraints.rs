use ark_poly::{Evaluations, Radix2EvaluationDomain};
use ark_bw6_761::Fr;
use crate::domains::Domains;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_ff::One;

/// Register polynomials in evaluation form amplified to support degree 4n constraints
pub(crate) struct Registers<'a> {
    domains: &'a Domains,
    bitmask: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    // // public keys' coordinates
    // pks_x: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    // pks_y: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    // // aggregate public key rolling sum coordinates
    // apk_acc_x: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    // apk_acc_y: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
}

impl<'a> Registers<'a> {
    pub fn new(bitmask: Vec<Fr>, domains: &'a Domains) -> Self {
        Self {
            domains,
            bitmask: domains.amplify(bitmask),
        }
    }
}

pub(crate) struct Constraints {}

impl Constraints {
    pub fn compute_bitmask_booleanity_constraint_polynomial(registers: &Registers) -> DensePolynomial<Fr> {
        let b = &registers.bitmask;
        let mut one_minus_b = registers.domains.constant_4x(Fr::one());
        one_minus_b -= b;
        (b * &one_minus_b).interpolate()
    }

    pub fn evaluate_bitmask_booleanity_constraint(bitmask_at_zeta: Fr) -> Fr {
        bitmask_at_zeta * (Fr::one() - bitmask_at_zeta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use ark_std::rand::{Rng, rngs::StdRng};
    use ark_ff::{Zero, Fp384};
    use ark_poly::Polynomial;
    use ark_bls12_377::FqParameters;

    fn random_bitmask(rng: &mut StdRng, n: usize) -> Vec<Fp384<FqParameters>> {
        (0..n)
            .map(|_| rng.gen_bool(2.0 / 3.0))
            .map(|b| if b { Fr::one() } else { Fr::zero() })
            .collect::<Vec<_>>()
    }

    #[test]
    fn test_bitmask_booleanity_constraint() {
        let rng = &mut test_rng();
        let n = 64;
        let domains = Domains::new(n);


        let good_bitmask = random_bitmask(rng, n);
        let registers = Registers::new(good_bitmask, &domains);
        let constraint_poly =
            Constraints::compute_bitmask_booleanity_constraint_polynomial(&registers);
        assert_eq!(constraint_poly.degree(), 2 * (n - 1));
        assert!(domains.is_zero(&constraint_poly));


        let mut bad_bitmask = random_bitmask(rng, n);
        bad_bitmask[0] = Fr::rand(rng);
        let registers = Registers::new(bad_bitmask, &domains);
        let constraint_poly =
            Constraints::compute_bitmask_booleanity_constraint_polynomial(&registers);
        assert_eq!(constraint_poly.degree(), 2 * (n - 1));
        assert!(!domains.is_zero(&constraint_poly));
    }
}