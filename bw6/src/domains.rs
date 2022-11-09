use ark_bls12_377::Fq;
use ark_bw6_761::Fr;
use ark_ff::{Field, One, Zero};
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Evaluations, Radix2EvaluationDomain};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_std::convert::TryInto;

#[derive(Clone)]
pub struct Domains {
    //TODO: remove pub
    pub domain: Radix2EvaluationDomain<Fr>,
    pub domain2x: Radix2EvaluationDomain<Fr>,
    pub domain4x: Radix2EvaluationDomain<Fr>,

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
        let domain2x =
            Radix2EvaluationDomain::<Fr>::new(2 * domain_size).unwrap();
        let domain4x =
            Radix2EvaluationDomain::<Fr>::new(4 * domain_size).unwrap();

        let l_first = Self::first_lagrange_basis_polynomial(domain_size);
        let l_last = Self::last_lagrange_basis_polynomial(domain_size);
        let l_first_evals_over_4x = Self::_amplify(l_first, domain, domain4x);
        let l_last_evals_over_4x = Self::_amplify(l_last, domain, domain4x);

        Domains {
            domain,
            domain2x,
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
        poly.evaluate_over_domain_by_ref(self.domain4x)
    }

    pub fn amplify(&self, evals: Vec<Fr>) -> Evaluations<Fr, Radix2EvaluationDomain<Fr>> {
        Self::_amplify(evals, self.domain, self.domain4x)
    }

    /// Checks if the polynomial is identically zero over the smaller domain.
    pub fn is_zero(&self, poly: &DensePolynomial<Fr>) -> bool {
        poly.divide_by_vanishing_poly(self.domain).unwrap().1 == DensePolynomial::zero()
    }

    /// Divides by the vanishing polynomial of the smaller domain.
    pub fn compute_quotient(&self, poly: &DensePolynomial<Fr>) -> (DensePolynomial<Fq>, DensePolynomial<Fq>) {
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

    pub fn amplify_x2(&self, evals: Vec<Fr>) -> Evaluations<Fr, Radix2EvaluationDomain<Fr>> {
        let evals = Evaluations::from_vec_and_domain(evals, self.domain);
        let poly = evals.interpolate_by_ref();
        let evals = evals.evals;

        let omega_2x = self.domain2x.group_gen;
        let coset_poly = Self::coset_polynomial(&poly, omega_2x);
        let coset_evals = coset_poly.evaluate_over_domain_by_ref(self.domain);
        let evals2x = evals.into_iter().zip(coset_evals.evals)
            .flat_map(|(e, ce)| vec![e, ce])
            .collect();
        Evaluations::from_vec_and_domain(evals2x, self.domain2x)
    }

    pub fn amplify_x4(&self, evals: Vec<Fr>) -> Evaluations<Fr, Radix2EvaluationDomain<Fr>> {
        let evals = Evaluations::from_vec_and_domain(evals, self.domain);
        let poly = evals.interpolate_by_ref();
        let evals = evals.evals;

        let omega_4x = self.domain4x.group_gen;
        let coset_evals: [Vec<Fr>; 3] = (1..4)
            .map(|i| omega_4x.pow([i]))
            .map(|gi| Self::coset_polynomial(&poly, gi))
            .map(|p| p.evaluate_over_domain_by_ref(self.domain).evals)
            .collect::<Vec<_>>().try_into().unwrap();
        let evals_4x = evals.iter()
            .zip(&coset_evals[0])
            .zip(&coset_evals[1])
            .zip(&coset_evals[2])
            .flat_map(|(((g0, g1), g2), g3)| [g0, g1, g2, g3])
            .cloned()
            .collect();
        Evaluations::from_vec_and_domain(evals_4x, self.domain4x)
    }

    /// For a polynomial p returns a polynomial p' such that p'(H) = p(gH)
    fn coset_polynomial(poly: &DensePolynomial<Fr>, g: Fr) -> DensePolynomial<Fr> {
        let coset_coeffs = poly.coeffs.iter()
            .scan(Fr::one(), |pow, &coeff| {
                let coset_coeff = *pow * coeff;
                *pow = *pow * g;
                Some(coset_coeff)
            })
            .collect();
        DensePolynomial::from_coefficients_vec(coset_coeffs)
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{test_rng, UniformRand};

    use super::*;

    #[test]
    fn test_coset_amplify() {
        use ark_poly::DenseUVPolynomial;

        let rng = &mut test_rng();
        let n = 64;

        // Let H < G be a subgroup of index 2 (meaning |G| = 2|H|).
        // Then G = H \cup gH, where g is a generator of G.

        // Let |H| = n + 1. Observe that
        // evaluations of a degree n polynomial p(X) = a_0 + ... + a_n.X^n over a coset gH
        // are equal to the
        // evaluations of the polynomial p'(X) = a_0 + ... + (a_n.g^n).X^n over the subgroup H:
        // p(gH) = p'(H).

        // Thus p(G) can be computed either with a 2n-FFT,
        // or as p(G) = p(H \cup gH) = p(H) \cup p(gH) = p(H) \cup p'(H) with 2 n-FFTs.
        // In the case when p(H) is already known, the latter approach might be more efficient.

        let domain = Radix2EvaluationDomain::<Fr>::new(n).unwrap(); // H
        let domain2x = Radix2EvaluationDomain::<Fr>::new(2 * n).unwrap(); // G
        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>(); // p(H)
        let poly = Evaluations::from_vec_and_domain(evals.clone(), domain).interpolate(); // p
        let evals2x = poly.evaluate_over_domain_by_ref(domain2x); // p(G)

        let root2x = domain2x.group_gen; // g
        let coset_coeffs = poly.coeffs.iter()
            .scan(Fr::one(), |pow, &coeff| {
                let coset_coeff = *pow * coeff;
                *pow = *pow * root2x;
                Some(coset_coeff)
            })
            .collect();

        let coset_poly = DensePolynomial::from_coefficients_vec(coset_coeffs); // p'
        let coset_evals = coset_poly.evaluate_over_domain_by_ref(domain); // p'(H)

        let evals2x_2: Vec<_> = evals.into_iter().zip(coset_evals.evals)
            .flat_map(|(e, ce)| vec![e, ce])
            .collect(); // p(G)

        assert_eq!(evals2x.evals, evals2x_2);
    }

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
    fn test_amplify_2x() {
        let rng = &mut test_rng();
        let n = 64;

        let domains = Domains::new(n);

        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let poly = Evaluations::from_vec_and_domain(evals.clone(), domains.domain).interpolate();
        let evals2x = poly.evaluate_over_domain_by_ref(domains.domain2x);

        let evals2x_2 = domains.amplify_x2(evals);

        assert_eq!(evals2x, evals2x_2);
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