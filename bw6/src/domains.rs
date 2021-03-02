use ark_poly::{Radix2EvaluationDomain, Evaluations, EvaluationDomain};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_ff::{Zero, One};
use ark_bw6_761::Fr;
use ark_bls12_377::Fq;

pub struct Domains {
    //TODO: remove pub
    pub domain: Radix2EvaluationDomain<Fr>, // TODO: separate type?
    domain4x: Radix2EvaluationDomain<Fr>, // TODO: separate type?

    /// First Lagrange basis polynomial L_0 of degree n evaluated over the domain of size 4 * n; L_0(\omega^0) = 1
    pub l_first_evals_over_4x: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    /// Last  Lagrange basis polynomial L_{n-1} of degree n evaluated over the domain of size 4 * n; L_{n-1}(\omega^{n-1}}) = 1
    pub l_last_evals_over_4x: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    /// \omega, a primitive n-th root of unity. Multiplicative generator of the smaller domain.
    pub omega: Fr,
    /// \omega^{n-1}
    pub omega_inv: Fr,
    /// The smaller domain size.
    pub size: usize,
}

impl Domains {
    pub fn new(domain_size: usize) -> Self {
        let domain =
            Radix2EvaluationDomain::<Fr>::new(domain_size).unwrap();
        let domain4x =
            Radix2EvaluationDomain::<Fr>::new(4 * domain_size).unwrap();

        let l_first = Self::first_lagrange_basis_polynomial(domain_size);
        let l_last = Self::last_lagrange_basis_polynomial(domain_size);
        let l_first_evals_over_4x = Self::_amplify(l_first, domain, domain4x);
        let l_last_evals_over_4x = Self::_amplify(l_last, domain, domain4x);

        Domains {
            domain,
            domain4x,
            l_first_evals_over_4x,
            l_last_evals_over_4x,
            omega: domain.group_gen,
            omega_inv: domain.group_gen_inv,
            size: domain.size(),
        }
    }

    /// Interpolates the evaluations over the smaller domain,
    /// resulting in a degree < n polynomial.
    pub fn interpolate(&self, evals: Vec<Fr>) -> DensePolynomial<Fr> {
        // TODO: assert evals.len()
        Evaluations::from_vec_and_domain(evals, self.domain).interpolate()
    }

    /// Produces evaluations of the degree < n polynomial over the larger domain,
    /// resulting in a vec of evaluations of length 4n.
    pub fn amplify_polynomial(&self, poly: &DensePolynomial<Fr>) -> Evaluations<Fr, Radix2EvaluationDomain<Fr>> {
        // TODO: assert poly.degree()
        poly.evaluate_over_domain_by_ref( self.domain4x)
    }

    pub fn amplify(&self, evals: Vec<Fr>) -> Evaluations<Fr, Radix2EvaluationDomain<Fr>> {
        Self::_amplify(evals, self.domain, self.domain4x)
    }

    /// Checks if the polynomial is identically zero over the smaller domain.
    pub fn is_zero(&self, poly: &DensePolynomial<Fr>) -> bool {
        poly.divide_by_vanishing_poly(self.domain).unwrap().1 == DensePolynomial::zero()
    }

    /// Divides by the vanishing polynomial of the smaller domain.
    pub fn compute_quotient(&self, poly: &DensePolynomial<Fr> ) -> (DensePolynomial<Fq>, DensePolynomial<Fq>) {
        poly.divide_by_vanishing_poly(self.domain).unwrap() //TODO: arkworks never returns None
    }

    /// Degree n polynomial c * L_{n-1} evaluated over domain of size 4 * n.
    pub fn l_last_scaled_by(&self, c: Fr) -> Evaluations<Fr, Radix2EvaluationDomain<Fr>> {
        &self.constant_4x(c) * &self.l_last_evals_over_4x
    }

    pub fn constant_4x(&self, c: Fr) -> Evaluations<Fr, Radix2EvaluationDomain<Fr>> {
        // TODO: ConstantEvaluations to save memory
        let evals = vec![c; self.domain4x.size()];
        Evaluations::from_vec_and_domain(evals, self.domain4x)
    }

    /// Produces evaluations of a degree n polynomial in 4n points, given evaluations in n points.
    /// That allows arithmetic operations with degree n polynomials in evaluations form until the result extends degree 4n.
    // TODO: test
    // takes nlogn + 4nlog(4n) = nlogn + 4nlogn + 8n
    // TODO: can we do better?
    fn _amplify(evals: Vec<Fr>, domain: Radix2EvaluationDomain<Fr>, domain4x: Radix2EvaluationDomain<Fr>) -> Evaluations<Fr, Radix2EvaluationDomain<Fr>> {
        let poly = Evaluations::from_vec_and_domain(evals, domain).interpolate();
        let evals4x = poly.evaluate_over_domain(domain4x);
        evals4x
    }

    fn first_lagrange_basis_polynomial(domain_size: usize) -> Vec<Fr> {
        Self::li(0, domain_size)
    }

    fn last_lagrange_basis_polynomial(domain_size: usize) -> Vec<Fr> {
        Self::li(domain_size - 1, domain_size)
    }

    fn li(i: usize, domain_size: usize) -> Vec<Fr> {
        let mut li = vec![Fr::zero(); domain_size];
        li[i] = Fr::one();
        li
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn test_amplify() {
        let rng = &mut test_rng();
        let n = 64;

        let domains = Domains::new(n);

        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let poly = domains.interpolate(evals.clone());

        let evals4x_from_poly = domains.amplify_polynomial(&poly);
        let evals4x_from_vec = domains.amplify(evals);

        assert_eq!(evals4x_from_poly, evals4x_from_vec);
        assert_eq!(evals4x_from_poly.interpolate(), poly);
    }

    #[test]
    fn test_domains_l_last_scaled_by() {
        let rng = &mut test_rng();
        let n = 64;

        let c = Fr::rand(rng);

        let mut c_ln = vec![Fr::zero(); n];
        c_ln[n - 1] = c;

        let domains = Domains::new(n);

        assert_eq!(domains.l_last_scaled_by(c), domains.amplify(c_ln));
    }
}