use ark_poly::{Evaluations, Radix2EvaluationDomain, UVPolynomial};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_ff::{One, Zero, Field, Fp384};
use ark_bw6_761::Fr;
use ark_bls12_377::{G1Affine, FqParameters};

use crate::domains::Domains;
use crate::{Bitmask, point_in_g1_complement, utils};
use ark_ec::short_weierstrass_jacobian::GroupAffine;
use crate::utils::LagrangeEvaluations;

/// Register polynomials in evaluation form amplified to support degree 4n constraints
#[derive(Clone)] //TODO: remove
pub(crate) struct Registers<'a> {
    domains: &'a Domains,
    bitmask: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    // public keys' coordinates
    pks_x: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    pks_y: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    // aggregate public key rolling sum coordinates
    apk_acc_x: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    apk_acc_y: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    // TODO: not really registers
    apk_acc_x_shifted: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    apk_acc_y_shifted: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    apk_plus_h: G1Affine,
}

impl<'a> Registers<'a> {
    pub fn new(domains: &'a Domains,
               bitmask: &Bitmask,
               pks: Vec<G1Affine>,
    ) -> Self {
        let m = pks.len();
        let n = domains.size;
        assert!(m + 1 <= n);  // keyset_size + 1 <= domain_size (accounts for partial sums acc initial value)

        assert_eq!(bitmask.size(), m);
        let h = point_in_g1_complement();
        let mut acc = vec![h; m + 1];
        for (i, (b, p)) in bitmask.to_bits().iter().zip(pks.clone()).enumerate() {
            acc[i + 1] = if *b {
                acc[i] + p
            } else {
                acc[i]
            }
        }

        let (mut acc_x, mut acc_y): (Vec<Fr>, Vec<Fr>) = acc.iter()
            .map(|p| (p.x, p.y))
            .unzip();

        assert_eq!(acc_x.len(), m + 1);
        assert_eq!(acc_y.len(), m + 1);
        assert_eq!(GroupAffine::new(acc_x[0], acc_y[0], false), h);
        // assert_eq!(GroupAffine::new(acc_x[m], acc_y[m], false), apk.into_affine() + h);

        let mut b = bitmask.to_bits().iter()
            .map(|b| if *b { Fr::one() } else { Fr::zero() })
            .collect::<Vec<_>>();

        // Extend the computation to the whole domain
        b.resize_with(n, || Fr::zero());
        // So we don't care about pks, but
        let apk_plus_h_x = acc_x[m];
        let apk_plus_h_y = acc_y[m];
        acc_x.resize_with(n, || apk_plus_h_x);
        acc_y.resize_with(n, || apk_plus_h_y);



        let mut bitmask = bitmask.to_bits_as_field_elements();
        bitmask.resize(domains.size, Fr::zero());

        let pks = pks.iter()
            .map(|p| (p.x, p.y))
            .unzip();

        let mut apk_acc_x_shifted = acc_x.clone();
        let mut apk_acc_y_shifted = acc_y.clone();
        apk_acc_x_shifted.rotate_left(1);
        apk_acc_y_shifted.rotate_left(1);

        Self::new_unchecked(
            domains,
            bitmask,
            pks,
            (acc_x, acc_y),
            (apk_acc_x_shifted, apk_acc_y_shifted),
        )
    }

    fn new_unchecked(domains: &'a Domains,
                     bitmask: Vec<Fr>,
                     pks: (Vec<Fr>, Vec<Fr>),
                     apk_acc: (Vec<Fr>, Vec<Fr>),
                     apk_acc_shifted: (Vec<Fr>, Vec<Fr>),
    ) -> Self {
        // TODO: lol this is really ugly
        let apk_plus_h = GroupAffine::new(*apk_acc.0.last().unwrap(), *apk_acc.1.last().unwrap(), false);
        Self {
            domains,
            bitmask: domains.amplify(bitmask),
            pks_x: domains.amplify(pks.0),
            pks_y: domains.amplify(pks.1),
            apk_acc_x: domains.amplify(apk_acc.0),
            apk_acc_y: domains.amplify(apk_acc.1),
            apk_acc_x_shifted: domains.amplify(apk_acc_shifted.0),
            apk_acc_y_shifted: domains.amplify(apk_acc_shifted.1),
            apk_plus_h,
        }
    }

    // TODO: interpolate over the smaller domain
    pub fn get_partial_sums_register_polynomials(&self) -> (DensePolynomial<Fr>, DensePolynomial<Fr>) {
        (self.apk_acc_x.interpolate_by_ref(), self.apk_acc_y.interpolate_by_ref())
    }

    // TODO: interpolate over the smaller domain
    pub fn get_bitmask_register_polynomial(&self) -> DensePolynomial<Fr> {
        self.bitmask.interpolate_by_ref()
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

    pub fn compute_conditional_affine_addition_constraint_polynomials(registers: &Registers) ->
    (DensePolynomial<Fr>, DensePolynomial<Fr>) {
        let b = &registers.bitmask;
        let mut one_minus_b = registers.domains.constant_4x(Fr::one());
        one_minus_b -= b;

        let x1 = &registers.apk_acc_x;
        let y1 = &registers.apk_acc_y;
        let x2 = &registers.pks_x;
        let y2 = &registers.pks_y;
        let x3 = &registers.apk_acc_x_shifted;
        let y3 = &registers.apk_acc_y_shifted;

        let c1 =
            &(
                b *
                    &(
                        &(
                            &(
                                &(x1 - x2) * &(x1 - x2)
                            ) *
                                &(
                                    &(x1 + x2) + x3
                                )
                        ) -
                            &(
                                &(y2 - y1) * &(y2 - y1)
                            )
                    )
            ) +
                &(
                    &one_minus_b * &(y3 - y1)
                );

        let c2 =
            &(
                b *
                    &(
                        &(
                            &(x1 - x2) * &(y3 + y1)
                        ) -
                            &(
                                &(y2 - y1) * &(x3 - x1)
                            )
                    )
            ) +
                &(
                    &one_minus_b * &(x3 - x1)
                );

        let c1_poly = c1.interpolate();
        let c2_poly = c2.interpolate();

        // Multiply by selector polynomial
        // ci *= (X - \omega^{n-1})
        let mut a1_poly_ = mul_by_x(&c1_poly);
        a1_poly_ += (-registers.domains.omega_inv, &c1_poly);
        let mut a2_poly_ = mul_by_x(&c2_poly);
        a2_poly_ += (-registers.domains.omega_inv, &c2_poly);
        (a1_poly_, a2_poly_)
    }

    pub fn evaluate_conditional_affine_addition_constraints(
        b: Fr,
        x1: Fr,
        y1: Fr,
        x2: Fr,
        y2: Fr,
        x3: Fr,
        y3: Fr,
    ) -> (Fr, Fr) {
        let c1 =
            b * (
                (x1 - x2) * (x1 - x2) * (x1 + x2 + x3)
                    - (y2 - y1) * (y2 - y1)
            ) + (Fr::one() - b) * (y3 - y1);

        let c2 =
            b * (
                (x1 - x2) * (y3 + y1)
                    - (y2 - y1) * (x3 - x1)
            ) + (Fr::one() - b) * (x3 - x1);

        (c1, c2)
    }

    pub fn evaluate_conditional_affine_addition_constraints_linearized(
        b: Fr,
        x1: Fr,
        y1: Fr,
        x2: Fr,
        y2: Fr,
    ) -> (Fr, Fr) {
        Self::evaluate_conditional_affine_addition_constraints(b, x1, y1, x2, y2, Fr::zero(), Fr::zero())
    }

    // TODO: better name
    pub fn compute_public_inputs_constraint_polynomials(registers: &Registers) ->
    (DensePolynomial<Fr>, DensePolynomial<Fr>) {
        let h = point_in_g1_complement();
        let x1 = &registers.apk_acc_x;
        let y1 = &registers.apk_acc_y;

        let acc_minus_h_x = x1 - &registers.domains.constant_4x(h.x);
        let acc_minus_h_y = y1 - &registers.domains.constant_4x(h.y);

        let acc_minus_h_plus_apk_x =
            x1 - &registers.domains.constant_4x(registers.apk_plus_h.x);
        let acc_minus_h_plus_apk_y =
            y1 - &registers.domains.constant_4x(registers.apk_plus_h.y);

        let a4 = &(&acc_minus_h_x * &registers.domains.l_first_evals_over_4x)
            + &(&acc_minus_h_plus_apk_x * &registers.domains.l_last_evals_over_4x);
        let a5 = &(&acc_minus_h_y * &registers.domains.l_first_evals_over_4x)
            + &(&acc_minus_h_plus_apk_y * &registers.domains.l_last_evals_over_4x);
        let a4_poly = a4.interpolate();
        let a5_poly = a5.interpolate();
        (a4_poly, a5_poly)
    }

    pub fn evaluate_public_inputs_constraints(
        apk: G1Affine,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
        x1: Fr,
        y1: Fr,
    ) -> (Fr, Fr) {
        let h = point_in_g1_complement();
        let apk_plus_h = h + apk;
        let c1 = (x1 - h.x) * evals_at_zeta.l_first + (x1 - apk_plus_h.x) * evals_at_zeta.l_last;
        let c2 = (y1 - h.y) * evals_at_zeta.l_first + (y1 - apk_plus_h.y) * evals_at_zeta.l_last;
        (c1, c2)
    }
}

pub(crate) struct SuccinctlyAccountableRegisters<'a> {
    registers: Registers<'a>,
    c: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    c_shifted: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    acc: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    acc_shifted: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
}

impl<'a> SuccinctlyAccountableRegisters<'a> {
    pub fn new(registers: Registers<'a>,
               bitmask: Vec<Fr>,
               bitmask_chunks_aggregation_challenge: Fr, // denoted 'r' in the write-ups
    ) -> Self {
        let n = registers.domains.size;
        let bits_in_bitmask_chunk = 256;  //256 is the highest power of 2 that fits field bit capacity //TODO: const
        assert_eq!(n % bits_in_bitmask_chunk, 0); // n is a power of 2

        let r = bitmask_chunks_aggregation_challenge;
        let c = Self::build_multipacking_mask_register(n, bits_in_bitmask_chunk, r);
        let acc = Self::build_partial_inner_products_register(n, &bitmask, &c);

        let mut c_shifted = c.clone();
        c_shifted.rotate_left(1);
        let mut acc_shifted = acc.clone();
        acc_shifted.rotate_left(1);

        Self::new_unchecked(
            registers,
            c,
            c_shifted,
            acc,
            acc_shifted,
        )
    }

    fn new_unchecked(
        registers: Registers<'a>,
        c: Vec<Fr>,
        c_shifted: Vec<Fr>,
        acc: Vec<Fr>,
        acc_shifted: Vec<Fr>,
    ) -> Self {
        Self {
            registers: registers.clone(), //TODO: fix
            c: registers.domains.amplify(c),
            c_shifted: registers.domains.amplify(c_shifted),
            acc: registers.domains.amplify(acc),
            acc_shifted: registers.domains.amplify(acc_shifted),
        }
    }

    //TODO: comment
    fn build_multipacking_mask_register(domain_size: usize, chunk_size: usize, randomizer: Fr) -> Vec<Fr> {
        let powers_of_2 = utils::powers(Fr::from(2u8), chunk_size - 1);
        let powers_of_r = utils::powers(randomizer, domain_size / chunk_size - 1);
        // tensor product (powers_of_r X powers_of_2)
        powers_of_r.iter().flat_map(|rj|
            powers_of_2.iter().map(move |_2k| *rj * _2k)
        ).collect::<Vec<Fr>>()
    }

    /// Returns length n vec (0, a[0]b[0],...,a[n-2]b[n-2]), where n is domain size
    fn build_partial_inner_products_register(domain_size: usize, a: &Vec<Fr>, b: &Vec<Fr>) -> Vec<Fr> {
        // we ignore the last elements but still...
        assert_eq!(a.len(), domain_size);
        assert_eq!(b.len(), domain_size);
        let mut acc = Vec::with_capacity(domain_size);
        acc.push(Fr::zero());
        a.iter().zip(b.iter())
            .map(|(a, b)| *a * b)
            .take(domain_size - 1)
            .for_each(|x| {
                acc.push(x + acc.last().unwrap());
            });
        acc
    }

    // TODO: interpolate over the smaller domain
    pub fn get_multipacking_mask_register_polynomial(&self) -> DensePolynomial<Fr> {
        self.c.interpolate_by_ref()
    }

    // TODO: interpolate over the smaller domain
    pub fn get_partial_inner_products_register_polynomial(&self) -> DensePolynomial<Fr> {
        self.acc.interpolate_by_ref()
    }
}

// TODO: implement multiplication by a sparse polynomial in arkworks?
fn mul_by_x<F: Field>(p: &DensePolynomial<F>) -> DensePolynomial<F> {
    let mut px = vec![F::zero()];
    px.extend_from_slice(&p.coeffs);
    DensePolynomial::from_coefficients_vec(px)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use ark_std::rand::{Rng, rngs::StdRng};
    use ark_poly::Polynomial;
    use ark_bls12_377::G1Projective;
    use ark_ec::{ProjectiveCurve, AffineCurve};
    use crate::tests::random_bits;
    use crate::utils;

    // TODO: there's crate::tests::random_bits
    fn random_bitmask(rng: &mut StdRng, n: usize) -> Vec<Fr> {
        (0..n)
            .map(|_| rng.gen_bool(2.0 / 3.0))
            .map(|b| if b { Fr::one() } else { Fr::zero() })
            .collect()
    }

    fn random_pks(n: usize, rng: &mut StdRng) -> Vec<G1Affine> {
        (0..n)
            .map(|_| G1Projective::rand(rng))
            .map(|p| p.into_affine())
            .collect()
    }

    fn dummy_registers(n: usize) -> (Vec<Fr>, Vec<Fr>) {
        (vec![Fr::zero(); n], vec![Fr::zero(); n])
    }

    #[test]
    fn test_bitmask_booleanity_constraint() {
        let rng = &mut test_rng();
        let n = 64;
        let m = n - 1;
        let domains = Domains::new(n);

        let good_bitmask = Bitmask::from_bits(&random_bits(m, 0.5, rng));
        let registers = Registers::new(
            &domains,
            &good_bitmask,
            random_pks(m, rng),
        );
        let constraint_poly =
            Constraints::compute_bitmask_booleanity_constraint_polynomial(&registers);
        assert_eq!(constraint_poly.degree(), 2 * (n - 1));
        assert!(domains.is_zero(&constraint_poly));
        let zeta = Fr::rand(rng);
        let bitmask_at_zeta = registers.get_bitmask_register_polynomial().evaluate(&zeta);
        assert_eq!(
            Constraints::evaluate_bitmask_booleanity_constraint(bitmask_at_zeta),
            constraint_poly.evaluate(&zeta)
        );

        let mut bad_bitmask = random_bitmask(rng, m);
        bad_bitmask[0] = Fr::rand(rng);
        let registers = Registers::new_unchecked(
            &domains,
            bad_bitmask,
            dummy_registers(n),
            dummy_registers(n),
            dummy_registers(n)
        );
        let constraint_poly =
            Constraints::compute_bitmask_booleanity_constraint_polynomial(&registers);
        assert_eq!(constraint_poly.degree(), 2 * (n - 1));
        assert!(!domains.is_zero(&constraint_poly));
    }

    #[test]
    fn test_conditional_affine_addition_constraints() {
        let rng = &mut test_rng();
        let n = 64;
        let m = n - 1;
        let domains = Domains::new(n);

        let bitmask = Bitmask::from_bits(&random_bits(m, 0.5, rng));
        let registers = Registers::new(
            &domains,
            &bitmask,
            random_pks(m, rng),
        );
        let constraint_polys =
            Constraints::compute_conditional_affine_addition_constraint_polynomials(&registers);
        assert_eq!(constraint_polys.0.degree(), 4 * n - 3);
        assert_eq!(constraint_polys.1.degree(), 3 * n - 2);
        assert!(domains.is_zero(&constraint_polys.0));
        assert!(domains.is_zero(&constraint_polys.1));
        // TODO: eval test
        // TODO: negative test?
    }

    #[test]
    fn test_public_inputs_constraints() {
        let rng = &mut test_rng();
        let n = 64;
        let m = n - 1;
        let domains = Domains::new(n);

        let bits = random_bits(m, 0.5, rng);
        let bitmask = Bitmask::from_bits(&bits);
        let pks = random_pks(m, rng);
        let registers = Registers::new(
            &domains,
            &bitmask,
            pks.clone(),
        );
        let constraint_polys =
            Constraints::compute_public_inputs_constraint_polynomials(&registers);
        assert_eq!(constraint_polys.0.degree(), 2 * n - 2);
        assert_eq!(constraint_polys.1.degree(), 2 * n - 2);
        assert!(domains.is_zero(&constraint_polys.0));
        assert!(domains.is_zero(&constraint_polys.1));

        //TODO: reuse smth else
        let apk = bits.iter()
            .zip(pks)
            .filter(|(&b, _p)| b)
            .map(|(_b, p)| p.into_projective())
            .sum::<ark_bls12_377::G1Projective>().into_affine();
        let zeta = Fr::rand(rng);
        let evals_at_zeta = utils::lagrange_evaluations(zeta, registers.domains.domain);
        let acc_polys = registers.get_partial_sums_register_polynomials();
        let (x1, y1) = (acc_polys.0.evaluate(&zeta), acc_polys.1.evaluate(&zeta));
        assert_eq!(
            Constraints::evaluate_public_inputs_constraints(apk, &evals_at_zeta, x1, y1),
            (constraint_polys.0.evaluate(&zeta), constraint_polys.1.evaluate(&zeta))
        );

        // TODO: negative test?
    }

    #[test]
    fn test_multipacking_mask_register() {
        let r = Fr::rand(&mut test_rng());
        let two = Fr::from(2u8);
        let multipacking_mask = SuccinctlyAccountableRegisters::build_multipacking_mask_register(4, 2, r);
        assert_eq!(multipacking_mask, vec![Fr::one(), two, r, r * two]);
    }

    #[test]
    fn test_partial_inner_products_register() {
        let from_u8_vec = |v: [u8; 4]| v.iter().map(|&x| Fr::from(x)).collect::<Vec<Fr>>();
        let a = from_u8_vec([1, 2, 3, 4]);
        let b = from_u8_vec([5, 6, 7, 8]);
        let partial_inner_product = SuccinctlyAccountableRegisters::build_partial_inner_products_register(4, &a, &b);
        assert_eq!(partial_inner_product, from_u8_vec([0, 1 * 5, 1 * 5 + 2 * 6, 1 * 5 + 2 * 6 + 3 * 7]));
    }


}