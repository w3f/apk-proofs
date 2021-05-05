use ark_poly::univariate::DensePolynomial;
use ark_bw6_761::{Fr, G1Affine};
use crate::piop::{RegisterPolynomials, RegisterEvaluations, RegisterPolys, RegisterCommitments};
use ark_poly::{Polynomial, Evaluations, Radix2EvaluationDomain, UVPolynomial};
use ark_ff::{Zero, One, Field};
use ark_ec::AffineCurve;
use crate::utils::LagrangeEvaluations;
use crate::{Bitmask, point_in_g1_complement};
use crate::domains::Domains;
use ark_ec::short_weierstrass_jacobian::GroupAffine;

use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct PartialSumsCommitments(pub ark_bw6_761::G1Affine, pub ark_bw6_761::G1Affine);

impl RegisterCommitments for PartialSumsCommitments {
    fn as_vec(&self) -> Vec<G1Affine> {
        vec![
            self.0,
            self.1,
        ]
    }
}

pub struct PartialSumsPolynomials(pub DensePolynomial<Fr>, pub DensePolynomial<Fr>);

impl RegisterPolys for PartialSumsPolynomials {
    type C = PartialSumsCommitments;

    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> PartialSumsCommitments {
        PartialSumsCommitments(f(&self.0), f(&self.1))
    }
}

#[derive(Clone)] //TODO: remove
pub struct BasicRegisterPolynomials {
    keyset: (DensePolynomial<Fr>, DensePolynomial<Fr>),
    bitmask: DensePolynomial<Fr>,
    partial_sums: (DensePolynomial<Fr>, DensePolynomial<Fr>),
}

impl RegisterPolynomials<AffineAdditionEvaluations> for BasicRegisterPolynomials {
    fn to_vec(self) -> Vec<DensePolynomial<Fr>> {
        vec![
            self.keyset.0,
            self.keyset.1,
            // self.bitmask,
            self.partial_sums.0,
            self.partial_sums.1,
        ]
    }

    fn evaluate(&self, point: Fr) -> AffineAdditionEvaluations {
        AffineAdditionEvaluations {
            keyset: (self.keyset.0.evaluate(&point), self.keyset.1.evaluate(&point)),
            bitmask: self.bitmask.evaluate(&point),
            partial_sums: (self.partial_sums.0.evaluate(&point), self.partial_sums.1.evaluate(&point)),
        }
    }
}


#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct AffineAdditionEvaluations {
    pub keyset: (Fr, Fr),
    pub bitmask: Fr,
    pub partial_sums: (Fr, Fr),
}

impl RegisterEvaluations for AffineAdditionEvaluations {
    type AC = ();
    type C = PartialSumsCommitments;

    fn as_vec(&self) -> Vec<Fr> {
        vec![
            self.keyset.0,
            self.keyset.1,
            // self.bitmask,
            self.partial_sums.0,
            self.partial_sums.1,
        ]
    }

    fn get_bitmask(&self) -> Fr {
        self.bitmask
    }

    fn restore_commitment_to_linearization_polynomial(&self,
                                                      phi: Fr,
                                                      zeta_minus_omega_inv: Fr,
                                                      commitments: &PartialSumsCommitments,
                                                      _extra_commitments: &(),
    ) -> ark_bw6_761::G1Projective {
        let b = self.bitmask;
        let (x1, y1) = self.partial_sums;
        let (x2, y2) = self.keyset;

        let mut r_comm = ark_bw6_761::G1Projective::zero();
        // X3 := acc_x polynomial
        // Y3 := acc_y polynomial
        // a1_lin = b(x1-x2)^2.X3 + (1-b)Y3
        // a2_lin = b(x1-x2)Y3 + b(y1-y2)X3 + (1-b)X3 // *= phi
        // X3 term = b(x1-x2)^2 + b(y1-y2)phi + (1-b)phi
        // Y3 term = (1-b) + b(x1-x2)phi
        // ...and both multiplied by (\zeta - \omega^{n-1}) // = zeta_minus_omega_inv
        r_comm += commitments.0.mul(zeta_minus_omega_inv * (b * (x1 - x2) * (x1 - x2) + b * (y1 - y2) * phi + (Fr::one() - b) * phi));
        r_comm += commitments.1.mul(zeta_minus_omega_inv * ((Fr::one() - b) + b * (x1 - x2) * phi));
        r_comm
    }

    fn evaluate_constraint_polynomials(
        &self,
        apk: ark_bls12_377::G1Affine,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
        r: Fr,
        bitmask: &Bitmask,
        domain_size: u64,
    ) -> Vec<Fr> {
        let b = self.bitmask;
        let (x1, y1) = self.partial_sums;
        let (x2, y2) = self.keyset;

        let (a1, a2) = Constraints::evaluate_conditional_affine_addition_constraints_linearized(evals_at_zeta.zeta_minus_omega_inv, b, x1, y1, x2, y2);
        let a3 = Constraints::evaluate_bitmask_booleanity_constraint(b);
        let (a4, a5) = Constraints::evaluate_public_inputs_constraints(apk, &evals_at_zeta, x1, y1);
        vec![a1, a2, a3, a4, a5]
    }

    fn is_accountable(&self) -> bool {
        false
    }
}

/// Register polynomials in evaluation form amplified to support degree 4n constraints
#[derive(Clone)] //TODO: remove
pub struct AffineAdditionRegisters {
    pub domains: Domains,
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
    apk_plus_h: ark_bls12_377::G1Affine,
    pub polynomials: BasicRegisterPolynomials,
}

impl AffineAdditionRegisters {
    pub fn new(domains: Domains,
               bitmask: &Bitmask,
               pks: Vec<ark_bls12_377::G1Affine>,
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

    fn new_unchecked(domains: Domains,
                     bitmask: Vec<Fr>,
                     pks: (Vec<Fr>, Vec<Fr>),
                     apk_acc: (Vec<Fr>, Vec<Fr>),
                     apk_acc_shifted: (Vec<Fr>, Vec<Fr>),
    ) -> Self {
        // TODO: lol this is really ugly
        let apk_plus_h = GroupAffine::new(*apk_acc.0.last().unwrap(), *apk_acc.1.last().unwrap(), false);

        let bitmask_polynomial = domains.interpolate(bitmask);
        let keyset_polynomial = (
            domains.interpolate(pks.0),
            domains.interpolate(pks.1),
        );
        let partial_sums_polynomial = (
            domains.interpolate(apk_acc.0),
            domains.interpolate(apk_acc.1),
        );

        Self {
            domains: domains.clone(),
            bitmask: domains.amplify_polynomial(&bitmask_polynomial),
            pks_x: domains.amplify_polynomial(&keyset_polynomial.0),
            pks_y: domains.amplify_polynomial(&keyset_polynomial.1),
            apk_acc_x: domains.amplify_polynomial(&partial_sums_polynomial.0),
            apk_acc_y: domains.amplify_polynomial(&partial_sums_polynomial.1),
            apk_acc_x_shifted: domains.amplify(apk_acc_shifted.0),
            apk_acc_y_shifted: domains.amplify(apk_acc_shifted.1),
            apk_plus_h,
            polynomials: BasicRegisterPolynomials {
                bitmask: bitmask_polynomial,
                keyset: keyset_polynomial,
                partial_sums: partial_sums_polynomial,
            },
        }
    }

    // TODO: interpolate over the smaller domain
    pub fn get_bitmask_register_polynomial(&self) -> DensePolynomial<Fr> {
        self.bitmask.interpolate_by_ref()
    }

    // TODO: interpolate over the smaller domain
    pub fn get_partial_sums_register_polynomials(&self) -> PartialSumsPolynomials {
        PartialSumsPolynomials(
            self.apk_acc_x.interpolate_by_ref(),
            self.apk_acc_y.interpolate_by_ref()
        )
    }

    pub fn evaluate_register_polynomials(&self, point: Fr) -> AffineAdditionEvaluations {
        self.polynomials.evaluate(point)
    }

    // Compute linearization polynomial
    // See https://hackmd.io/CdZkCe2PQuy7XG7CLOBRbA step 4
    // deg(r) = n, so it can be computed in the monomial basis
    pub fn compute_linearization_polynomial(&self, evaluations: &AffineAdditionEvaluations, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr> {
        let b_zeta = evaluations.bitmask;
        let (acc_x_zeta, acc_y_zeta) = (evaluations.partial_sums.0, evaluations.partial_sums.1);
        let (pks_x_zeta, pks_y_zeta) = (evaluations.keyset.0, evaluations.keyset.1);
        let (acc_x_poly, acc_y_poly) = &self.polynomials.partial_sums;

        let mut a1_lin = DensePolynomial::<Fr>::zero();
        a1_lin += (b_zeta * (acc_x_zeta - pks_x_zeta) * (acc_x_zeta - pks_x_zeta), acc_x_poly);
        a1_lin += (Fr::one() - b_zeta, acc_y_poly);

        let mut a2_lin = DensePolynomial::<Fr>::zero();
        a2_lin += (b_zeta * (acc_x_zeta - pks_x_zeta), acc_y_poly);
        a2_lin += (b_zeta * (acc_y_zeta - pks_y_zeta), acc_x_poly);
        a2_lin += (Fr::one() - b_zeta, acc_x_poly);

        let mut r_poly = DensePolynomial::<Fr>::zero();
        r_poly += (zeta_minus_omega_inv, &a1_lin);
        r_poly += (zeta_minus_omega_inv * phi, &a2_lin);
        r_poly
    }

    pub fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>> {
        let (a1_poly, a2_poly) =
            Constraints::compute_conditional_affine_addition_constraint_polynomials(self);
        let a3_poly =
            Constraints::compute_bitmask_booleanity_constraint_polynomial(self);
        let (a4_poly, a5_poly) =
            Constraints::compute_public_inputs_constraint_polynomials(self);
        vec![a1_poly, a2_poly, a3_poly, a4_poly, a5_poly]
    }

    pub fn get_all_register_polynomials(self) -> Vec<DensePolynomial<Fr>> {
        self.polynomials.to_vec()
    }
}

pub(crate) struct Constraints {}

impl Constraints {
    pub fn compute_bitmask_booleanity_constraint_polynomial(registers: &AffineAdditionRegisters) -> DensePolynomial<Fr> {
        let b = &registers.bitmask;
        let mut one_minus_b = registers.domains.constant_4x(Fr::one());
        one_minus_b -= b;
        (b * &one_minus_b).interpolate()
    }

    pub fn evaluate_bitmask_booleanity_constraint(bitmask_at_zeta: Fr) -> Fr {
        bitmask_at_zeta * (Fr::one() - bitmask_at_zeta)
    }

    pub fn compute_conditional_affine_addition_constraint_polynomials(registers: &AffineAdditionRegisters) ->
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
        zeta_minus_omega_inv: Fr,
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

        (c1 * zeta_minus_omega_inv, c2 * zeta_minus_omega_inv)
    }

    pub fn evaluate_conditional_affine_addition_constraints_linearized(
        zeta_minus_omega_inv: Fr,
        b: Fr,
        x1: Fr,
        y1: Fr,
        x2: Fr,
        y2: Fr,
    ) -> (Fr, Fr) {
        Self::evaluate_conditional_affine_addition_constraints(zeta_minus_omega_inv, b, x1, y1, x2, y2, Fr::zero(), Fr::zero())
    }

    // TODO: better name
    pub fn compute_public_inputs_constraint_polynomials(registers: &AffineAdditionRegisters) ->
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
        apk: ark_bls12_377::G1Affine,
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

    fn random_pks(n: usize, rng: &mut StdRng) -> Vec<ark_bls12_377::G1Affine> {
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
        let registers = AffineAdditionRegisters::new(
            domains.clone(),
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
        let registers = AffineAdditionRegisters::new_unchecked(
            domains.clone(),
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
        let registers = AffineAdditionRegisters::new(
            domains.clone(),
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
        let registers = AffineAdditionRegisters::new(
            domains.clone(),
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
}

