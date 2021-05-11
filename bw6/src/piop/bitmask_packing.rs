use ark_poly::{Evaluations, Radix2EvaluationDomain, Polynomial};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_ff::{One, Zero, Field};
use ark_bw6_761::Fr;
use ark_ec::AffineCurve;

use ark_std::io::{Read, Write};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};
use ark_std::{end_timer, start_timer};

use crate::{Bitmask, utils};
use crate::utils::LagrangeEvaluations;
use crate::piop::{VerifierProtocol, RegisterCommitments, RegisterPolynomials, RegisterEvaluations};
use crate::piop::affine_addition::{AffineAdditionEvaluations, PartialSumsAndBitmaskCommitments};
use crate::domains::Domains;


#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct BitmaskPackingCommitments {
    pub c_comm: ark_bw6_761::G1Affine,
    pub acc_comm: ark_bw6_761::G1Affine,
}

impl BitmaskPackingCommitments {
    pub fn new(c_comm: ark_bw6_761::G1Affine, acc_comm: ark_bw6_761::G1Affine) -> Self {
        BitmaskPackingCommitments { c_comm, acc_comm }
    }
}

impl RegisterCommitments for BitmaskPackingCommitments {
    fn as_vec(&self) -> Vec<ark_bw6_761::G1Affine> {
        vec![
            self.c_comm,
            self.acc_comm,
        ]
    }
}

#[derive(Clone)]
pub struct BitmaskPackingPolynomials {
    pub c_poly: DensePolynomial<Fr>,
    pub acc_poly: DensePolynomial<Fr>,
}

impl BitmaskPackingPolynomials {
    pub fn new(c_poly: DensePolynomial<Fr>, acc_poly: DensePolynomial<Fr>) -> Self {
        BitmaskPackingPolynomials { c_poly, acc_poly }
    }

    //TODO: &self
    pub fn to_vec(self) -> Vec<DensePolynomial<Fr>> {
        vec![
            self.c_poly,
            self.acc_poly,
        ]
    }
}

impl RegisterPolynomials for BitmaskPackingPolynomials {
    type C = BitmaskPackingCommitments;

    fn commit<F: Fn(&DensePolynomial<Fr>) -> ark_bw6_761::G1Affine>(&self, f: F) -> Self::C {
        BitmaskPackingCommitments::new(
            f(&self.c_poly),
            f(&self.acc_poly),
        )
    }
}

//TODO: remove pubs
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct SuccinctAccountableRegisterEvaluations {
    pub c: Fr,
    pub acc: Fr,
    pub basic_evaluations: AffineAdditionEvaluations,
}

impl RegisterEvaluations for SuccinctAccountableRegisterEvaluations {
    fn as_vec(&self) -> Vec<Fr> {
        let mut res = self.basic_evaluations.as_vec();
        res.extend(vec![self.c, self.acc]);
        res
    }
}

impl VerifierProtocol for SuccinctAccountableRegisterEvaluations {
    type AC = BitmaskPackingCommitments;
    type C = PartialSumsAndBitmaskCommitments;



    fn restore_commitment_to_linearization_polynomial(&self,
                                                      phi: Fr,
                                                      zeta_minus_omega_inv: Fr,
                                                      commitments: &PartialSumsAndBitmaskCommitments,
                                                      extra_commitments: &BitmaskPackingCommitments,
    ) -> ark_bw6_761::G1Projective {
        let powers_of_phi = utils::powers(phi, 6);
        let mut r_comm = self.basic_evaluations.restore_commitment_to_linearization_polynomial(phi, zeta_minus_omega_inv, &commitments.partial_sums, &());
        r_comm += extra_commitments.acc_comm.mul(powers_of_phi[5]);
        r_comm += extra_commitments.c_comm.mul(powers_of_phi[6]);
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

        let a6 = BitmaskPackingRegisters::evaluate_inner_product_constraint_linearized(
            aggregated_bitmask,
            &evals_at_zeta,
            b,
            c,
            acc,
        );

        let a7 = BitmaskPackingRegisters::evaluate_multipacking_mask_constraint_linearized(
            a,
            r_pow_m,
            &evals_at_zeta,
            c,
        );

        let mut res = self.basic_evaluations.evaluate_constraint_polynomials(apk, evals_at_zeta, r, bitmask, domain_size);
        res.extend(vec![a6, a7]);
        res
    }
}



pub(crate) struct BitmaskPackingRegisters {
    domains: Domains,

    bitmask: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    c: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    c_shifted: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    acc: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    acc_shifted: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,

    bitmask_chunks_aggregated: Fr,
    polynomials: BitmaskPackingPolynomials,
    r: Fr,
}

impl BitmaskPackingRegisters {

    // TODO: remove bitmask arg
    pub fn new(domains: Domains,
               bitmask: &Bitmask,
               bitmask_chunks_aggregation_challenge: Fr, // denoted 'r' in the write-ups
    ) -> Self {
        let n = domains.size;
        let bits_in_bitmask_chunk = 256;  //256 is the highest power of 2 that fits field bit capacity //TODO: const
        assert_eq!(n % bits_in_bitmask_chunk, 0); // n is a power of 2

        let mut bitmask = bitmask.to_bits_as_field_elements();
        bitmask.resize(domains.size, Fr::zero());

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
            domains,
            bitmask,
            c,
            c_shifted,
            acc,
            acc_shifted,
            bitmask_chunks_aggregated,
            r
        )
    }

    fn new_unchecked(
        domains: Domains,

        bitmask: Vec<Fr>,
        c: Vec<Fr>,
        c_shifted: Vec<Fr>,
        acc: Vec<Fr>,
        acc_shifted: Vec<Fr>,
        bitmask_chunks_aggregated: Fr,
        r: Fr
    ) -> Self {
        let c_polynomial = domains.interpolate(c);
        let acc_polynomial = domains.interpolate(acc);
        Self {
            domains: domains.clone(),

            bitmask: domains.amplify(bitmask),
            c: domains.amplify_polynomial(&c_polynomial),
            c_shifted: domains.amplify(c_shifted),
            acc: domains.amplify_polynomial(&acc_polynomial),
            acc_shifted: domains.amplify(acc_shifted),
            bitmask_chunks_aggregated,
            polynomials: BitmaskPackingPolynomials {
                c_poly: c_polynomial,
                acc_poly: acc_polynomial,
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
        let bc_ln_x4 = self.domains.l_last_scaled_by(self.bitmask_chunks_aggregated);
        let constraint = &(&(&self.acc_shifted - &self.acc) - &(&self.bitmask * &self.c)) + &bc_ln_x4;
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
        let n = self.domains.size;
        let chunks = n / 256; //TODO: consts
        let mut a = vec![Fr::from(2u8); n];
        a.iter_mut().step_by(256).for_each(|a| *a = self.r / Fr::from(2u8).pow([255u64]));
        a.rotate_left(1);
        let a_x4 = self.domains.amplify(a);

        let x_todo = Fr::one() - self.r.pow([chunks as u64]); //TODO: name
        let ln_x4 = self.domains.l_last_scaled_by(x_todo);

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



impl BitmaskPackingRegisters {
    pub fn evaluate_register_polynomials(&self, point: Fr) -> (Fr, Fr) {
        //TODO: struct
        (
            self.polynomials.c_poly.evaluate(&point),
            self.polynomials.acc_poly.evaluate(&point),
        )
    }

    pub fn compute_constraints_linearized(&self) -> Vec<DensePolynomial<Fr>> {
        vec![
            self.polynomials.acc_poly.clone(),
            self.polynomials.c_poly.clone(),
        ]
    }

    pub fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>> {
        let a6_poly = self.compute_inner_product_constraint_polynomial();
        let a7_poly = self.compute_multipacking_mask_constraint_polynomial();
        vec![a6_poly, a7_poly]
    }

    pub fn get_register_polynomials(&self) -> BitmaskPackingPolynomials {
        self.polynomials.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use ark_std::rand::{Rng, rngs::StdRng};
    use ark_poly::Polynomial;
    use ark_bls12_377::G1Projective;
    use ark_ec::ProjectiveCurve;
    use crate::tests::random_bits;
    use crate::domains::Domains;
    use crate::piop::affine_addition::AffineAdditionRegisters;

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
    fn test_multipacking_mask_register() {
        let r = Fr::rand(&mut test_rng());
        let two = Fr::from(2u8);
        let multipacking_mask = BitmaskPackingRegisters::build_multipacking_mask_register(4, 2, r);
        assert_eq!(multipacking_mask, vec![Fr::one(), two, r, r * two]);
    }

    #[test]
    fn test_partial_inner_products_register() {
        let from_u8_vec = |v: [u8; 4]| v.iter().map(|&x| Fr::from(x)).collect::<Vec<Fr>>();
        let a = from_u8_vec([1, 2, 3, 4]);
        let b = from_u8_vec([5, 6, 7, 8]);
        let partial_inner_product = BitmaskPackingRegisters::build_partial_inner_products_register(4, &a, &b);
        assert_eq!(partial_inner_product, from_u8_vec([0, 1 * 5, 1 * 5 + 2 * 6, 1 * 5 + 2 * 6 + 3 * 7]));
    }

    #[test]
    fn test_multipacking_mask_constraint() {
        let rng = &mut test_rng();
        let n = 256;
        let m = n - 1;
        let domains = Domains::new(n);

        let bitmask = Bitmask::from_bits(&random_bits(m, 0.5, rng));

        let r = Fr::rand(rng);
        let acc_registers = BitmaskPackingRegisters::new(
            domains.clone(),
            &bitmask,
            r
        );
        let constraint_poly = acc_registers.compute_multipacking_mask_constraint_polynomial();
        assert_eq!(constraint_poly.degree(), 2 * n - 2);
        assert!(domains.is_zero(&constraint_poly));
    }
}