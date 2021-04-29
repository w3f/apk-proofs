use ark_poly::{Evaluations, Radix2EvaluationDomain, UVPolynomial, Polynomial};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_ff::{One, Zero, Field};
use ark_bw6_761::Fr;
use ark_ec::AffineCurve;
use ark_ec::short_weierstrass_jacobian::GroupAffine;
use ark_bls12_377::{G1Affine, Fq};

use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};
use ark_std::{end_timer, start_timer};

use crate::domains::Domains;
use crate::{Bitmask, point_in_g1_complement, utils, RegisterCommitments};
use crate::utils::LagrangeEvaluations;
use crate::piop::{Protocol, RegisterPolynomials, PackedAccountabilityRegisterPolynomials, PartialSumsPolynomials, PartialSumsCommitments, PackedRegisterCommitments, RegisterPolys};

#[derive(Clone)] //TODO: remove
pub struct BasicRegisterPolynomials {
    keyset: (DensePolynomial<Fr>, DensePolynomial<Fr>),
    bitmask: DensePolynomial<Fr>,
    partial_sums: (DensePolynomial<Fr>, DensePolynomial<Fr>),
}

impl RegisterPolynomials<BasicRegisterEvaluations> for BasicRegisterPolynomials {
    fn to_vec(self) -> Vec<DensePolynomial<Fr>> {
        vec![
            self.keyset.0,
            self.keyset.1,
            // self.bitmask,
            self.partial_sums.0,
            self.partial_sums.1,
        ]
    }

    fn evaluate(&self, point: Fr) -> BasicRegisterEvaluations {
        BasicRegisterEvaluations {
            keyset: (self.keyset.0.evaluate(&point), self.keyset.1.evaluate(&point)),
            bitmask: self.bitmask.evaluate(&point),
            partial_sums: (self.partial_sums.0.evaluate(&point), self.partial_sums.1.evaluate(&point)),
        }
    }
}

pub(crate) struct SuccinctAccountableRegisterPolynomials {
    basic_polynomials: BasicRegisterPolynomials,
    pub c_poly: DensePolynomial<Fr>,
    pub acc_poly: DensePolynomial<Fr>,
}

impl RegisterPolynomials<SuccinctAccountableRegisterEvaluations> for SuccinctAccountableRegisterPolynomials {
    fn to_vec(self) -> Vec<DensePolynomial<Fr>> {
        let mut res = self.basic_polynomials.to_vec();
        res.extend(vec![self.c_poly, self.acc_poly]);
        res
    }

    fn evaluate(&self, point: Fr) -> SuccinctAccountableRegisterEvaluations {
        SuccinctAccountableRegisterEvaluations {
            c: self.c_poly.evaluate(&point),
            acc: self.acc_poly.evaluate(&point),
            basic_evaluations: self.basic_polynomials.evaluate(point)
        }
    }
}

pub trait RegisterEvaluations {
    type AC: RegisterCommitments;
    type C: RegisterCommitments;

    fn as_vec(&self) -> Vec<Fr>;
    fn get_bitmask(&self) -> Fr;
    fn restore_commitment_to_linearization_polynomial(
        &self,
        phi: Fr,
        zeta_minus_omega_inv: Fr,
        commitments: &Self::C,
        extra_commitments: &Self::AC,
    ) -> ark_bw6_761::G1Projective;

    fn evaluate_constraint_polynomials(
        &self,
        apk: G1Affine,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
        r: Fr,
        bitmask: &Bitmask,
        domain_size: u64,
    ) -> Vec<Fr>;

    //TODO: move somewhere
    fn is_accountable(&self) -> bool;
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct BasicRegisterEvaluations {
    pub keyset: (Fr, Fr),
    pub bitmask: Fr,
    pub partial_sums: (Fr, Fr),
}

impl RegisterEvaluations for BasicRegisterEvaluations {
    type AC = ();
    type C = PartialSumsCommitments;

    fn as_vec(&self) -> Vec<Fq> {
        vec![
            self.keyset.0,
            self.keyset.1,
            // self.bitmask,
            self.partial_sums.0,
            self.partial_sums.1,
        ]
    }

    fn get_bitmask(&self) -> Fq {
        self.bitmask
    }

    fn restore_commitment_to_linearization_polynomial(&self,
                                                          phi: Fr,
                                                          zeta_minus_omega_inv: Fr,
                                                          commitments: &PartialSumsCommitments,
                                                          extra_commitments: &(),
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
        apk: G1Affine,
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

//TODO: remove pubs
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct SuccinctAccountableRegisterEvaluations {
    pub c: Fr,
    pub acc: Fr,
    pub basic_evaluations: BasicRegisterEvaluations,
}

impl RegisterEvaluations for SuccinctAccountableRegisterEvaluations {
    type AC = PackedRegisterCommitments;
    type C = PartialSumsCommitments;

    fn as_vec(&self) -> Vec<Fq> {
        let mut res = self.basic_evaluations.as_vec();
        res.extend(vec![self.c, self.acc]);
        res
    }

    fn get_bitmask(&self) -> Fq {
        self.basic_evaluations.bitmask
    }

    fn restore_commitment_to_linearization_polynomial(&self,
                                                      phi: Fr,
                                                      zeta_minus_omega_inv: Fr,
                                                      commitments: &PartialSumsCommitments,
                                                      extra_commitments: &PackedRegisterCommitments,
    ) -> ark_bw6_761::G1Projective {
        let powers_of_phi = utils::powers(phi, 6);
        let mut r_comm = self.basic_evaluations.restore_commitment_to_linearization_polynomial(phi, zeta_minus_omega_inv, &commitments, &());
        r_comm += extra_commitments.acc_comm.mul(powers_of_phi[5]);
        r_comm += extra_commitments.c_comm.mul(powers_of_phi[6]);
        r_comm
    }

    fn evaluate_constraint_polynomials(
        &self,
        apk: G1Affine,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
        r: Fr,
        bitmask: &Bitmask,
        domain_size: u64,
    ) -> Vec<Fr> {
        let bits_in_bitmask_chunk = 256;
        let bits_in_big_int_limb = 64;
        assert_eq!(bits_in_bitmask_chunk % bits_in_big_int_limb, 0);
        let limbs_in_chunk = bits_in_bitmask_chunk / bits_in_big_int_limb;
        assert_eq!(domain_size % bits_in_bitmask_chunk, 0);
        let chunks_in_bitmask = domain_size / bits_in_bitmask_chunk; // TODO: bitmask should be right-padded with 0s to domain_size

        let bits_in_bitmask_chunk_inv = Fr::from(256u16).inverse().unwrap();

        let powers_of_r = utils::powers(r, (chunks_in_bitmask - 1) as usize);
        let r_pow_m = r * powers_of_r.last().unwrap();
        let mut bitmask_chunks = bitmask.to_chunks_as_field_elements::<Fr>(limbs_in_chunk as usize);
        //TODO: pad in Bitmask
        bitmask_chunks.resize_with(chunks_in_bitmask as usize, || Fr::zero());
        assert_eq!(powers_of_r.len(), bitmask_chunks.len());
        let aggregated_bitmask = bitmask_chunks.into_iter()
            .zip(powers_of_r)
            .map(|(bj, rj)| bj * rj)
            .sum::<Fr>();


        let t_a_zeta_omega1 = start_timer!(|| "A(zw) as fraction");
        let zeta_omega_pow_m = evals_at_zeta.zeta_omega.pow([chunks_in_bitmask]); // m = chunks_in_bitmask
        let zeta_omega_pow_n = zeta_omega_pow_m.pow([bits_in_bitmask_chunk]); // n = domain_size
        let a_zeta_omega1 = bits_in_bitmask_chunk_inv * (zeta_omega_pow_n - Fr::one()) / (zeta_omega_pow_m - Fr::one());
        end_timer!(t_a_zeta_omega1);

        let t_a_zeta_omega2 = start_timer!(|| "A(zw) as polynomial");
        let zeta_omega_pow_m = evals_at_zeta.zeta_omega.pow([chunks_in_bitmask]); // m = chunks_in_bitmask
        let a_zeta_omega2 = bits_in_bitmask_chunk_inv * utils::powers(zeta_omega_pow_m, (bits_in_bitmask_chunk - 1) as usize).iter().sum::<Fr>();
        end_timer!(t_a_zeta_omega2);

        assert_eq!(a_zeta_omega1, a_zeta_omega2);
        let two = Fr::from(2u8);
        let a = two + (r / two.pow([255u64]) - two) * a_zeta_omega1;


        let b = self.basic_evaluations.bitmask;
        let acc = self.acc;
        let c = self.c;

        let a6 = SuccinctlyAccountableRegisters::evaluate_inner_product_constraint_linearized(
            aggregated_bitmask,
            &evals_at_zeta,
            b,
            c,
            acc,
        );

        let a7 = SuccinctlyAccountableRegisters::evaluate_multipacking_mask_constraint_linearized(
            a,
            r_pow_m,
            &evals_at_zeta,
            c,
        );

        let mut res = self.basic_evaluations.evaluate_constraint_polynomials(apk, evals_at_zeta, r, bitmask, domain_size);
        res.extend(vec![a6, a7]);
        res
    }

    fn is_accountable(&self) -> bool {
        true
    }
}


/// Register polynomials in evaluation form amplified to support degree 4n constraints
#[derive(Clone)] //TODO: remove
pub struct Registers {
    domains: Domains,
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
    pub polynomials: BasicRegisterPolynomials,
}

impl Registers {
    pub fn new(domains: Domains,
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

pub(crate) struct SuccinctlyAccountableRegisters {
    registers: Registers,
    c: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    c_shifted: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    acc: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    acc_shifted: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,

    bitmask_chunks_aggregated: Fr,
    pub polynomials: SuccinctAccountableRegisterPolynomials,
    r: Fr,
}

impl SuccinctlyAccountableRegisters {

    // TODO: remove bitmask arg
    pub fn new(registers: Registers,
               bitmask: Vec<Fr>,
               bitmask_chunks_aggregation_challenge: Fr, // denoted 'r' in the write-ups
    ) -> Self {
        let n = registers.domains.size;
        let bits_in_bitmask_chunk = 256;  //256 is the highest power of 2 that fits field bit capacity //TODO: const
        assert_eq!(n % bits_in_bitmask_chunk, 0); // n is a power of 2

        let r = bitmask_chunks_aggregation_challenge;
        let c = Self::build_multipacking_mask_register(n, bits_in_bitmask_chunk, r);
        let acc = Self::build_partial_inner_products_register(n, &bitmask, &c);
        let bitmask_chunks_aggregated = bitmask.iter()
            .zip(c.iter())
            .map(|(&b, c)| b * c)
            .sum::<Fr>();

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
            bitmask_chunks_aggregated,
            r
        )
    }

    fn new_unchecked(
        registers: Registers,
        c: Vec<Fr>,
        c_shifted: Vec<Fr>,
        acc: Vec<Fr>,
        acc_shifted: Vec<Fr>,
        bitmask_chunks_aggregated: Fr,
        r: Fr
    ) -> Self {
        let c_polynomial = registers.domains.interpolate(c);
        let acc_polynomial = registers.domains.interpolate(acc);
        Self {
            registers: registers.clone(), //TODO: fix
            c: registers.domains.amplify_polynomial(&c_polynomial),
            c_shifted: registers.domains.amplify(c_shifted),
            acc: registers.domains.amplify_polynomial(&acc_polynomial),
            acc_shifted: registers.domains.amplify(acc_shifted),
            bitmask_chunks_aggregated,
            polynomials: SuccinctAccountableRegisterPolynomials {
                c_poly: c_polynomial,
                acc_poly: acc_polynomial,
                basic_polynomials: registers.polynomials
            },
            r
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

    pub fn compute_inner_product_constraint_polynomial(&self) -> DensePolynomial<Fr> {
        let bc_ln_x4 = self.registers.domains.l_last_scaled_by(self.bitmask_chunks_aggregated);
        let constraint = &(&(&self.acc_shifted - &self.acc) - &(&self.registers.bitmask * &self.c)) + &bc_ln_x4;
        constraint.interpolate()
    }

    pub fn evaluate_inner_product_constraint(
        bitmask_chunks_aggregated: Fr,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
        b_zeta: Fr,
        c_zeta: Fr,
        acc_zeta: Fr,
        acc_zeta_omega: Fr,
    ) -> Fr {
        acc_zeta_omega - acc_zeta - b_zeta * c_zeta + bitmask_chunks_aggregated * evals_at_zeta.l_last
    }

    pub fn evaluate_inner_product_constraint_linearized(
        bitmask_chunks_aggregated: Fr,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
        b_zeta: Fr,
        c_zeta: Fr,
        acc_zeta: Fr
    ) -> Fr {
        Self::evaluate_inner_product_constraint(bitmask_chunks_aggregated, evals_at_zeta, b_zeta, c_zeta, acc_zeta, Fr::zero())
    }

    pub fn compute_multipacking_mask_constraint_polynomial(&self) -> DensePolynomial<Fr> {
        let n = self.registers.domains.size;
        let chunks = n / 256; //TODO: consts
        let mut a = vec![Fr::from(2u8); n];
        a.iter_mut().step_by(256).for_each(|a| *a = self.r / Fr::from(2u8).pow([255u64]));
        a.rotate_left(1);
        let a_x4 = self.registers.domains.amplify(a);

        let x_todo = Fr::one() - self.r.pow([chunks as u64]); //TODO: name
        let ln_x4 = self.registers.domains.l_last_scaled_by(x_todo);

        let a7 = &(&self.c_shifted - &(&self.c * &a_x4)) - &ln_x4;
        a7.interpolate()
    }

    pub fn evaluate_multipacking_mask_constraint(
        a: Fr,
        r_pow_m: Fr,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
        c_zeta: Fr,
        c_zeta_omega: Fr
    ) -> Fr {
        c_zeta_omega - c_zeta * a - (Fr::one() - r_pow_m) * evals_at_zeta.l_last
    }

    pub fn evaluate_multipacking_mask_constraint_linearized(
        a: Fr,
        r_pow_m: Fr,
        evals_at_zeta: &LagrangeEvaluations<Fr>,
        c_zeta: Fr,
    ) -> Fr {
        Self::evaluate_multipacking_mask_constraint(a, r_pow_m, evals_at_zeta, c_zeta, Fr::zero())
    }
}



impl  SuccinctlyAccountableRegisters {
    pub fn evaluate_register_polynomials(&self, point: Fr) -> SuccinctAccountableRegisterEvaluations {
        self.polynomials.evaluate(point)
    }

    pub fn compute_linearization_polynomial(&self, evaluations: &SuccinctAccountableRegisterEvaluations, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr> {
        let powers_of_phi = &utils::powers(phi, 6);
        // let a6 = &(&(&acc_shifted_x4 - &acc_x4) - &(&B * &c_x4)) + &(bc_ln_x4);
        let a6_lin = &self.polynomials.acc_poly;
        // let a7 = &(&c_shifted_x4 - &(&c_x4 * &a_x4)) - &ln_x4;
        let a7_lin = &self.polynomials.c_poly;

        let mut r_poly = self.registers.compute_linearization_polynomial(&evaluations.basic_evaluations, phi, zeta_minus_omega_inv);
        r_poly += (powers_of_phi[5], a6_lin);
        r_poly += (powers_of_phi[6], a7_lin);
        r_poly
    }

    pub fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fq>> {
        let mut constraints = self.registers.compute_constraint_polynomials();
        let a6_poly = self.compute_inner_product_constraint_polynomial();
        let a7_poly = self.compute_multipacking_mask_constraint_polynomial();
        constraints.extend(vec![a6_poly, a7_poly]);
        constraints
    }

    pub fn get_all_register_polynomials(self) -> Vec<DensePolynomial<Fr>> {
        self.polynomials.to_vec()
    }
}

impl Protocol<BasicRegisterEvaluations> for Registers {
    type P1 = PartialSumsPolynomials;
    type P2 = ();

    fn init(domains: Domains, bitmask: &Bitmask, pks: Vec<G1Affine>) -> Self {
        Registers::new(domains, bitmask, pks)
    }

    // TODO: interpolate over the smaller domain
    fn get_1st_round_register_polynomials(&self) -> Self::P1 {
        PartialSumsPolynomials(
            self.apk_acc_x.interpolate_by_ref(),
            self.apk_acc_y.interpolate_by_ref()
        )
    }

    fn evaluate_register_polynomials(&self, point: Fr) -> BasicRegisterEvaluations {
        self.polynomials.evaluate(point)
    }

    // Compute linearization polynomial
    // See https://hackmd.io/CdZkCe2PQuy7XG7CLOBRbA step 4
    // deg(r) = n, so it can be computed in the monomial basis
    fn compute_linearization_polynomial(&self, evaluations: &BasicRegisterEvaluations, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr> {
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

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fq>> {
        let (a1_poly, a2_poly) =
            Constraints::compute_conditional_affine_addition_constraint_polynomials(self);
        let a3_poly =
            Constraints::compute_bitmask_booleanity_constraint_polynomial(self);
        let (a4_poly, a5_poly) =
            Constraints::compute_public_inputs_constraint_polynomials(self);
        vec![a1_poly, a2_poly, a3_poly, a4_poly, a5_poly]
    }

    fn get_all_register_polynomials(self) -> Vec<DensePolynomial<Fq>> {
        self.polynomials.to_vec()
    }

    fn get_2nd_round_register_polynomials(&mut self, bitmask: Vec<Fr>, verifier_challenge: Fr) -> () {
        ()
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
        let registers = Registers::new_unchecked(
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
        let registers = Registers::new(
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
        let registers = Registers::new(
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
        let acc_polys = registers.get_1st_round_register_polynomials();
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

    #[test]
    fn test_multipacking_mask_constraint() {
        let rng = &mut test_rng();
        let n = 256;
        let m = n - 1;
        let domains = Domains::new(n);

        let bitmask = Bitmask::from_bits(&random_bits(m, 0.5, rng));
        let registers = Registers::new(
            domains.clone(),
            &bitmask,
            random_pks(m, rng),
        );

        let mut b = bitmask.to_bits_as_field_elements();
        b.resize_with(n, || Fr::zero());
        let r = Fr::rand(rng);
        let acc_registers = SuccinctlyAccountableRegisters::new(
            registers,
            b,
            r
        );
        let constraint_poly = acc_registers.compute_multipacking_mask_constraint_polynomial();
        assert_eq!(constraint_poly.degree(), 2 * n - 2);
        assert!(domains.is_zero(&constraint_poly));
    }
}