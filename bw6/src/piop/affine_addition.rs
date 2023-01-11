use std::iter;

use ark_bls12_377::G1Projective;
use ark_bw6_761::{Fr, G1Affine};
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{Field, One, Zero};
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Evaluations, Polynomial, Radix2EvaluationDomain};
use ark_poly::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use crate::{Keyset, point_in_g1_complement};
use crate::domains::Domains;
use crate::piop::{RegisterCommitments, RegisterEvaluations, RegisterPolynomials, VerifierProtocol};
use crate::utils::LagrangeEvaluations;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct PartialSumsCommitments(
    pub ark_bw6_761::G1Affine,
    pub ark_bw6_761::G1Affine,
);

impl RegisterCommitments for PartialSumsCommitments {
    fn as_vec(&self) -> Vec<G1Affine> {
        vec![
            self.0,
            self.1,
        ]
    }
}

pub type PartialSumsPolynomials = [DensePolynomial<Fr>; 2];

impl RegisterPolynomials for PartialSumsPolynomials {
    type C = PartialSumsCommitments;

    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> PartialSumsCommitments {
        PartialSumsCommitments(f(&self[0]), f(&self[1]))
    }
}

//TODO: move to packed.rs?

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct PartialSumsAndBitmaskCommitments {
    pub partial_sums: PartialSumsCommitments,
    pub bitmask: ark_bw6_761::G1Affine,
}

impl RegisterCommitments for PartialSumsAndBitmaskCommitments {
    fn as_vec(&self) -> Vec<G1Affine> {
        let mut res = vec![self.bitmask];
        res.extend(self.partial_sums.as_vec());
        res
    }
}

pub struct PartialSumsAndBitmaskPolynomials {
    pub partial_sums: PartialSumsPolynomials,
    pub bitmask: DensePolynomial<Fr>,
}

impl RegisterPolynomials for PartialSumsAndBitmaskPolynomials {
    type C = PartialSumsAndBitmaskCommitments;

    fn commit<F: Clone + Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> PartialSumsAndBitmaskCommitments {
        PartialSumsAndBitmaskCommitments {
            partial_sums: self.partial_sums.commit(f.clone()),
            bitmask: f(&self.bitmask),
        }
    }
}

#[derive(Clone)] //TODO: remove
pub struct AffineAdditionPolynomials {
    pub keyset: [DensePolynomial<Fr>; 2],
    pub bitmask: DensePolynomial<Fr>,
    pub partial_sums: [DensePolynomial<Fr>; 2],
}

impl AffineAdditionPolynomials {
    pub fn to_vec(self) -> Vec<DensePolynomial<Fr>> {
        IntoIterator::into_iter(self.keyset)
            .chain(std::iter::once(self.bitmask))
            .chain(self.partial_sums)
            .collect()
    }

    fn evaluate(&self, point: Fr) -> AffineAdditionEvaluations {
        AffineAdditionEvaluations {
            keyset: (self.keyset[0].evaluate(&point), self.keyset[1].evaluate(&point)),
            bitmask: self.bitmask.evaluate(&point),
            partial_sums: (self.partial_sums[0].evaluate(&point), self.partial_sums[1].evaluate(&point)),
        }
    }
}


#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct AffineAdditionEvaluations {
    pub keyset: (Fr, Fr),
    pub bitmask: Fr,
    pub partial_sums: (Fr, Fr),
}

impl RegisterEvaluations for AffineAdditionEvaluations {
    fn as_vec(&self) -> Vec<Fr> {
        vec![
            self.keyset.0,
            self.keyset.1,
            self.bitmask,
            self.partial_sums.0,
            self.partial_sums.1,
        ]
    }
}

impl VerifierProtocol for AffineAdditionEvaluations {
    type C2 = ();
    type C1 = PartialSumsCommitments;

    const POLYS_OPENED_AT_ZETA: usize = 5;

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
        r_comm += commitments.0 * (zeta_minus_omega_inv * (b * (x1 - x2) * (x1 - x2) + b * (y1 - y2) * phi + (Fr::one() - b) * phi));
        r_comm += commitments.1 * (zeta_minus_omega_inv * ((Fr::one() - b) + b * (x1 - x2) * phi));
        r_comm
    }
}

impl AffineAdditionEvaluations {
    pub fn evaluate_constraint_polynomials(
        &self,
        apk: ark_bls12_377::G1Affine,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
    ) -> Vec<Fr> {
        let b = self.bitmask;
        let (x1, y1) = self.partial_sums;
        let (x2, y2) = self.keyset;

        let (a1, a2) = Constraints::evaluate_conditional_affine_addition_constraints_linearized(evals_at_zeta.zeta_minus_omega_inv, b, x1, y1, x2, y2);
        let a3 = Constraints::evaluate_bitmask_booleanity_constraint(b);
        let (a4, a5) = Constraints::evaluate_public_inputs_constraints(apk, &evals_at_zeta, x1, y1);
        vec![a1, a2, a3, a4, a5]
    }
}

/// Register polynomials in evaluation form amplified to support degree 4n constraints
pub struct AffineAdditionRegisters {
    pub domains: Domains,
    bitmask: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    // public keys' coordinates
    keyset: [Evaluations<Fr, Radix2EvaluationDomain<Fr>>; 2],
    // aggregate public key rolling sum coordinates
    partial_sums: [Evaluations<Fr, Radix2EvaluationDomain<Fr>>; 2],

    pub polynomials: AffineAdditionPolynomials,
}

impl AffineAdditionRegisters {
    pub fn new(domains: Domains,
               keyset: Keyset,
               bitmask: &[bool],
    ) -> Self {
        assert_eq!(bitmask.len(), keyset.size());
        let domain_size = keyset.domain.size();

        let h = point_in_g1_complement().into_group();
        let apk_acc = bitmask.iter().zip(keyset.pks.iter())
            .scan(h, |acc, (b, pk)| {
                if *b {
                    *acc += pk;
                }
                Some(*acc)
            });
        let apk_acc: Vec<_> = iter::once(h)
            .chain(apk_acc)
            .collect();
        let mut apk_acc = G1Projective::normalize_batch(&apk_acc);

        apk_acc.resize(domain_size, apk_acc.last().cloned().unwrap());
        let apk_acc = apk_acc.iter()
            .map(|p| (p.x, p.y))
            .unzip();

        let mut bitmask = bitmask.to_vec();
        bitmask.resize(domain_size - 1, false);

        let bitmask = bitmask.iter()
            .map(|b| if *b { Fr::one() } else { Fr::zero() })
            .chain(iter::once(Fr::zero())) //TODO: pad with Fr::one()
            .collect();

        Self::new_unchecked(
            domains,
            bitmask,
            keyset,
            [apk_acc.0, apk_acc.1],
        )
    }

    fn new_unchecked(domains: Domains,
                     bitmask: Vec<Fr>,
                     keyset: Keyset,
                     apk_acc: [Vec<Fr>; 2],
    ) -> Self {
        let bitmask_polynomial = domains.interpolate(bitmask);
        let partial_sums_polynomial = apk_acc.map(|z| domains.interpolate(z));
        let partial_sums = partial_sums_polynomial.clone().map(|z| domains.amplify_polynomial(&z));
        let bitmask = domains.amplify_polynomial(&bitmask_polynomial);

        Self {
            domains,
            bitmask,
            keyset: keyset.pks_evals_x4.unwrap(),
            partial_sums,
            polynomials: AffineAdditionPolynomials {
                bitmask: bitmask_polynomial,
                keyset: keyset.pks_polys,
                partial_sums: partial_sums_polynomial,
            },
        }
    }

    pub fn evaluate_register_polynomials(&self, point: Fr) -> AffineAdditionEvaluations {
        self.polynomials.evaluate(point)
    }

    // Compute linearization polynomial
    // See https://hackmd.io/CdZkCe2PQuy7XG7CLOBRbA step 4
    // deg(r) = n, so it can be computed in the monomial basis
    pub fn compute_constraints_linearized(&self, evaluations: &AffineAdditionEvaluations, zeta: Fr) -> Vec<DensePolynomial<Fr>> {
        let zeta_minus_omega_inv = zeta - self.domains.omega_inv;
        let b_zeta = evaluations.bitmask;
        let (acc_x_zeta, acc_y_zeta) = (evaluations.partial_sums.0, evaluations.partial_sums.1);
        let (pks_x_zeta, pks_y_zeta) = (evaluations.keyset.0, evaluations.keyset.1);
        let [acc_x_poly, acc_y_poly] = &self.polynomials.partial_sums;

        let mut a1_lin = DensePolynomial::<Fr>::zero();
        a1_lin += (b_zeta * (acc_x_zeta - pks_x_zeta) * (acc_x_zeta - pks_x_zeta), acc_x_poly);
        a1_lin += (Fr::one() - b_zeta, acc_y_poly);
        // a1_lin = zeta_minus_omega_inv * a1_lin // TODO: fix in arkworks
        a1_lin.coeffs.iter_mut().for_each(|c| *c *= zeta_minus_omega_inv);

        let mut a2_lin = DensePolynomial::<Fr>::zero();
        a2_lin += (b_zeta * (acc_x_zeta - pks_x_zeta), acc_y_poly);
        a2_lin += (b_zeta * (acc_y_zeta - pks_y_zeta), acc_x_poly);
        a2_lin += (Fr::one() - b_zeta, acc_x_poly);
        // a2_lin = zeta_minus_omega_inv * a2_lin // TODO: fix in arkworks
        a2_lin.coeffs.iter_mut().for_each(|c| *c *= zeta_minus_omega_inv);

        vec![
            a1_lin,
            a2_lin,
            DensePolynomial::zero(),
            DensePolynomial::zero(),
            DensePolynomial::zero(),
        ]
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

    pub fn get_register_polynomials(&self) -> AffineAdditionPolynomials {
        self.polynomials.clone()
    }

    pub fn get_partial_sums_and_bitmask_polynomials(&self) -> PartialSumsAndBitmaskPolynomials {
        let polys = self.get_register_polynomials();
        PartialSumsAndBitmaskPolynomials {
            partial_sums: polys.partial_sums,
            bitmask: polys.bitmask,
        }
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

        let [x1, y1] = &registers.partial_sums;
        let [x2, y2] = &registers.keyset;
        let mut next_partial_sums = registers.partial_sums.clone();
        next_partial_sums.iter_mut().for_each(|z| z.evals.rotate_left(4));
        let [x3, y3] = &next_partial_sums;

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
        let [x1, y1] = &registers.partial_sums;
        let [h_x, h_y] = [x1, y1].map(|z| z[0]);
        let [apk_plus_h_x, apk_plus_h_y] = [x1, y1].map(|z| z[4 * (registers.domains.size - 1)]);

        let acc_minus_h_x = x1 - &registers.domains.constant_4x(h_x);
        let acc_minus_h_y = y1 - &registers.domains.constant_4x(h_y);

        let acc_minus_h_plus_apk_x =
            x1 - &registers.domains.constant_4x(apk_plus_h_x);
        let acc_minus_h_plus_apk_y =
            y1 - &registers.domains.constant_4x(apk_plus_h_y);

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
        let apk_plus_h = (h + apk).into_affine();
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
    use ark_ec::CurveGroup;
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};

    use crate::test_helpers::{_random_bitmask, _random_bits, random_pks};
    use crate::utils;

    use super::*;

    fn dummy_registers(n: usize) -> [Vec<Fr>; 2] {
        [vec![Fr::one(); n], vec![Fr::one(); n]]
    }

    #[test]
    fn test_bitmask_booleanity_constraint() {
        let rng = &mut test_rng();
        let n = 64;
        let m = n - 1;
        let domains = Domains::new(n);

        let good_bitmask = _random_bits(m, 0.5, rng);
        let mut keyset = Keyset::new(random_pks(m, rng));
        keyset.amplify();
        let registers = AffineAdditionRegisters::new(
            domains.clone(),
            keyset.clone(),
            &good_bitmask,
        );
        let constraint_poly =
            Constraints::compute_bitmask_booleanity_constraint_polynomial(&registers);
        assert_eq!(constraint_poly.degree(), 2 * (n - 1));
        assert!(domains.is_zero(&constraint_poly));
        let zeta = Fr::rand(rng);
        let bitmask_at_zeta = registers.get_register_polynomials().bitmask.evaluate(&zeta);
        assert_eq!(
            Constraints::evaluate_bitmask_booleanity_constraint(bitmask_at_zeta),
            constraint_poly.evaluate(&zeta)
        );

        let mut bad_bitmask = _random_bitmask(m, rng);
        bad_bitmask[0] = Fr::rand(rng);

        let registers = AffineAdditionRegisters::new_unchecked(
            domains.clone(),
            bad_bitmask,
            keyset,
            dummy_registers(n),
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

        let mut keyset = Keyset::new(random_pks(m, rng));
        keyset.amplify();
        let registers = AffineAdditionRegisters::new(
            domains.clone(),
            keyset,
            &_random_bits(m, 0.5, rng),
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

        let bits = _random_bits(m, 0.5, rng);

        let mut keyset = Keyset::new(random_pks(m, rng));
        keyset.amplify();
        let registers = AffineAdditionRegisters::new(
            domains.clone(),
            keyset.clone(),
            &bits,
        );
        let constraint_polys =
            Constraints::compute_public_inputs_constraint_polynomials(&registers);
        assert_eq!(constraint_polys.0.degree(), 2 * n - 2);
        assert_eq!(constraint_polys.1.degree(), 2 * n - 2);
        assert!(domains.is_zero(&constraint_polys.0));
        assert!(domains.is_zero(&constraint_polys.1));

        let apk = keyset.aggregate(&bits).into_affine();
        let zeta = Fr::rand(rng);
        let evals_at_zeta = utils::lagrange_evaluations(zeta, registers.domains.domain);
        let acc_polys = registers.get_register_polynomials().partial_sums;
        let (x1, y1) = (acc_polys[0].evaluate(&zeta), acc_polys[1].evaluate(&zeta));
        assert_eq!(
            Constraints::evaluate_public_inputs_constraints(apk, &evals_at_zeta, x1, y1),
            (constraint_polys.0.evaluate(&zeta), constraint_polys.1.evaluate(&zeta))
        );

        // TODO: negative test?
    }
}

