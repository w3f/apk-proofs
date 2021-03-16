use ark_bw6_761::{BW6_761, Fr};
use ark_ec::ProjectiveCurve;
use ark_ec::short_weierstrass_jacobian::GroupAffine;
use ark_ff::{Field, One, Zero};
use ark_poly::{Evaluations, Polynomial, UVPolynomial, Radix2EvaluationDomain};
use ark_poly::univariate::DensePolynomial;
use merlin::Transcript;

use crate::{KZG_BW6, Proof, point_in_g1_complement, Bitmask, utils};
use crate::transcript::ApkTranscript;
use crate::signer_set::SignerSetCommitment;
use crate::kzg::ProverKey;
use crate::bls::PublicKey;
use crate::domains::Domains;
use crate::constraints::{Registers, Constraints, SuccinctlyAccountableRegisters};




struct Params {
    domain_size: usize,
    kzg_pk: ProverKey<BW6_761>,
    h: ark_bls12_377::G1Affine,
}


struct Session<'a> {
    pks: &'a [PublicKey],
    pks_x_poly: DensePolynomial<Fr>,
    pks_y_poly: DensePolynomial<Fr>,
    pks_x_poly_evals_x4: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
    pks_y_poly_evals_x4: Evaluations<Fr, Radix2EvaluationDomain<Fr>>,
}


impl<'a> Session<'a> {
    pub fn new(pks: &'a[PublicKey], domains: &Domains) -> Self {
        let (pks_x, pks_y): (Vec<Fr>, Vec<Fr>) = pks.iter()
            .map(|p| p.0.into_affine())
            .map(|p| (p.x, p.y))
            .unzip();
        let pks_x_poly = domains.interpolate(pks_x);
        let pks_y_poly = domains.interpolate(pks_y);
        let pks_x_poly_evals_x4 = domains.amplify_polynomial(&pks_x_poly);
        let pks_y_poly_evals_x4 = domains.amplify_polynomial(&pks_y_poly);
        Self {
            pks,
            pks_x_poly,
            pks_y_poly,
            pks_x_poly_evals_x4,
            pks_y_poly_evals_x4,
        }
    }

    fn compute_apk(&self, bitmask: &[bool]) -> ark_bls12_377::G1Projective {
        bitmask.iter()
            .zip(self.pks.iter())
            .filter(|(b, _p)| **b)
            .map(|(_b, p)| p.0)
            .sum::<ark_bls12_377::G1Projective>()
    }
}


pub struct Prover<'a> {
    params: Params,
    domains: Domains,
    session: Session<'a>,
    preprocessed_transcript: Transcript,
}


impl<'a> Prover<'a> {
    pub fn new(
        domain_size: usize,
        kzg_pk: ProverKey<BW6_761>,
        signer_set_comm: &SignerSetCommitment,
        pks: &'a [PublicKey],
        mut empty_transcript: Transcript,
    ) -> Self {
        assert!(domain_size.is_power_of_two(), "domain size should be a power of 2");
        assert!(domain_size <= kzg_pk.max_coeffs(), "domain size shouldn't exceed srs length");

        // empty_transcript.set_protocol_params(); //TODO
        empty_transcript.set_signer_set(&signer_set_comm);

        let params = Params {
            domain_size,
            kzg_pk,
            h: point_in_g1_complement(),
        };

        let domains = Domains::new(domain_size);
        let session = Session::new(pks, &domains);

        Self {
            params,
            domains,
            session,
            preprocessed_transcript: empty_transcript,
        }
    }

    #[allow(non_snake_case)]
    pub fn prove(&self, bitmask: &Bitmask) -> Proof {
        let m = self.session.pks.len();
        let n = self.params.domain_size;

        assert_eq!(bitmask.size(), m);
        assert!(bitmask.count_ones() > 0);

        let apk = self.session.compute_apk(&bitmask.to_bits());

        let mut transcript = self.preprocessed_transcript.clone();
        transcript.append_public_input(&apk.into(), bitmask);


        let mut b = bitmask.to_bits().iter()
            .map(|b| if *b { Fr::one() } else { Fr::zero() })
            .collect::<Vec<_>>();

        // Extend the computation to the whole domain
        b.resize_with(n, || Fr::zero());

        // TODO: move to Session
        let pks = self.session.pks.iter()
            .map(|p| p.0.into_affine())
            .collect();

        let registers = Registers::new(&self.domains, bitmask, pks);

        let b_poly = registers.get_bitmask_register_polynomial();
        let (acc_x_poly, acc_y_poly) = registers.get_partial_sums_register_polynomials();

        let b_comm = KZG_BW6::commit(&self.params.kzg_pk, &b_poly);
        let acc_x_comm = KZG_BW6::commit(&self.params.kzg_pk, &acc_x_poly);
        let acc_y_comm = KZG_BW6::commit(&self.params.kzg_pk, &acc_y_poly);

        assert_eq!(b_poly.coeffs.len(), n);
        assert_eq!(b_poly.degree(), n - 1);

        let (a1_poly, a2_poly) =
            Constraints::compute_conditional_affine_addition_constraint_polynomials(&registers);
        let a3_poly =
            Constraints::compute_bitmask_booleanity_constraint_polynomial(&registers);
        let (a4_poly, a5_poly) =
            Constraints::compute_public_inputs_constraint_polynomials(&registers);

        transcript.append_proof_point(b"b_comm", &b_comm);
        transcript.append_proof_point(b"acc_x_comm", &acc_x_comm);
        transcript.append_proof_point(b"acc_y_comm", &acc_y_comm);
        let r = transcript.get_128_bit_challenge(b"r"); // bitmask batching challenge

        let acc_registers = SuccinctlyAccountableRegisters::new(registers, b, r);
        let a6_poly = acc_registers.compute_inner_product_constraint_polynomial();
        let a7_poly = acc_registers.compute_multipacking_mask_constraint_polynomial(r);

        let c_poly = acc_registers.get_multipacking_mask_register_polynomial();
        let acc_poly = acc_registers.get_partial_inner_products_register_polynomial();

        let c_comm = KZG_BW6::commit(&self.params.kzg_pk, &c_poly);
        let acc_comm = KZG_BW6::commit(&self.params.kzg_pk, &acc_poly);
        transcript.append_proof_point(b"c_comm", &c_comm);
        transcript.append_proof_point(b"acc_comm", &acc_comm);
        let phi = transcript.get_128_bit_challenge(b"phi"); // constraint polynomials batching challenge

        let powers_of_phi = &utils::powers(phi, 6);
        let w = utils::randomize(phi, &[
            a1_poly,
            a2_poly,
            a3_poly,
            a4_poly,
            a5_poly,
            a6_poly,
            a7_poly,
        ]);

        let (q_poly, r) = self.domains.compute_quotient(&w);
        assert_eq!(r, DensePolynomial::zero());
        assert_eq!(q_poly.degree(), 3 * n - 3);

        assert_eq!(self.params.kzg_pk.max_degree(), q_poly.degree()); //TODO: check at the prover creation
        let q_comm = KZG_BW6::commit(&self.params.kzg_pk, &q_poly);

        transcript.append_proof_point(b"q_comm", &q_comm);
        let zeta = transcript.get_128_bit_challenge(b"zeta"); // evaluation point challenge

        let b_zeta = b_poly.evaluate(&zeta);
        let pks_x_zeta = self.session.pks_x_poly.evaluate(&zeta);
        let pks_y_zeta = self.session.pks_y_poly.evaluate(&zeta);
        let acc_x_zeta = acc_x_poly.evaluate(&zeta);
        let acc_y_zeta = acc_y_poly.evaluate(&zeta);
        let q_zeta = q_poly.evaluate(&zeta);
        let c_zeta = c_poly.evaluate(&zeta);
        let acc_zeta = acc_poly.evaluate(&zeta);

        let zeta_omega = zeta * self.domains.omega;
        let zeta_minus_omega_inv = zeta - self.domains.omega_inv;

        // Compute linearization polynomial
        // See https://hackmd.io/CdZkCe2PQuy7XG7CLOBRbA step 4
        // deg(r) = n, so it can be computed in the monomial basis

        let mut a1_lin = DensePolynomial::<Fr>::zero();
        a1_lin += (b_zeta * (acc_x_zeta - pks_x_zeta) * (acc_x_zeta - pks_x_zeta), &acc_x_poly);
        a1_lin += (Fr::one() - b_zeta, &acc_y_poly);

        let mut a2_lin = DensePolynomial::<Fr>::zero();
        a2_lin += (b_zeta * (acc_x_zeta - pks_x_zeta), &acc_y_poly);
        a2_lin += (b_zeta * (acc_y_zeta - pks_y_zeta), &acc_x_poly);
        a2_lin += (Fr::one() - b_zeta, &acc_x_poly);

        // let a6 = &(&(&acc_shifted_x4 - &acc_x4) - &(&B * &c_x4)) + &(bc_ln_x4);
        let a6_lin = acc_poly.clone();
        // let a7 = &(&c_shifted_x4 - &(&c_x4 * &a_x4)) - &ln_x4;
        let a7_lin = c_poly.clone();

        let mut r_poly = DensePolynomial::<Fr>::zero();
        r_poly += (zeta_minus_omega_inv, &a1_lin);
        r_poly += (zeta_minus_omega_inv * powers_of_phi[1], &a2_lin);
        r_poly += (powers_of_phi[5], &a6_lin);
        r_poly += (powers_of_phi[6], &a7_lin);

        let r_zeta_omega = r_poly.evaluate(&zeta_omega);


        transcript.append_proof_scalar(b"b_zeta", &b_zeta);
        transcript.append_proof_scalar(b"pks_x_zeta", &pks_x_zeta);
        transcript.append_proof_scalar(b"pks_y_zeta", &pks_y_zeta);
        transcript.append_proof_scalar(b"acc_x_zeta", &acc_x_zeta);
        transcript.append_proof_scalar(b"acc_y_zeta", &acc_y_zeta);
        transcript.append_proof_scalar(b"q_zeta", &q_zeta);
        transcript.append_proof_scalar(b"c_zeta", &c_zeta);
        transcript.append_proof_scalar(b"acc_zeta", &acc_zeta);
        transcript.append_proof_scalar(b"r_zeta_omega", &r_zeta_omega);
        let nu: Fr = transcript.get_128_bit_challenge(b"nu"); // KZG opening batching challenge

        let w_poly = KZG_BW6::aggregate_polynomials(nu, &[
            self.session.pks_x_poly.clone(),
            self.session.pks_y_poly.clone(),
            b_poly,
            q_poly,
            acc_poly,
            c_poly,
            acc_x_poly,
            acc_y_poly
        ]);

        let w_at_zeta_proof = KZG_BW6::open(&self.params.kzg_pk, &w_poly, zeta);
        let r_at_zeta_omega_proof = KZG_BW6::open(&self.params.kzg_pk, &r_poly, zeta_omega);

        transcript.append_proof_point(b"w_at_zeta_proof", &w_at_zeta_proof);
        transcript.append_proof_point(b"r_at_zeta_omega_proof", &r_at_zeta_omega_proof);

        Proof {
            b_comm,
            acc_x_comm,
            acc_y_comm,
            // r <-
            c_comm,
            acc_comm,
            // phi <-
            q_comm,
            // zeta <-
            b_zeta,
            pks_x_zeta,
            pks_y_zeta,
            acc_x_zeta,
            acc_y_zeta,
            q_zeta,
            c_zeta,
            acc_zeta,
            r_zeta_omega,
            // <- nu
            w_at_zeta_proof,
            r_at_zeta_omega_proof,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};


    #[test]
    fn test_mul_domain() {
        let rng = &mut test_rng();
        let n = 2;
        let d1 = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let d4 = GeneralEvaluationDomain::<Fr>::new(4 * n).unwrap();

        let a_evals1 = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let b_evals1 = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let a_poly1 = Evaluations::from_vec_and_domain(a_evals1, d1).interpolate();
        let b_poly1 = Evaluations::from_vec_and_domain(b_evals1, d1).interpolate();
        assert_eq!(a_poly1.degree(), n - 1);
        assert_eq!(b_poly1.degree(), n - 1);

        let a_evals4 = a_poly1.evaluate_over_domain_by_ref(d4);
        let b_evals4 = b_poly1.evaluate_over_domain_by_ref(d4);

        let c_evals4 = &a_evals4 * &b_evals4;
        let c_poly4 = c_evals4.interpolate();

        assert_eq!(c_poly4.degree(), 2 * (n - 1));
    }

    #[test]
    fn test_shift() {
        let rng = &mut test_rng();
        let n = 2;
        let d = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let p_evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let mut p_evals_shifted = p_evals.clone();
        p_evals_shifted.rotate_left(1);

        let p = Evaluations::from_vec_and_domain(p_evals, d).interpolate();
        let p_shifted = Evaluations::from_vec_and_domain(p_evals_shifted, d).interpolate();

        if let ark_poly::GeneralEvaluationDomain::Radix2(d) = d {
            let omega = d.group_gen;
            assert_eq!(p.evaluate(&omega), p_shifted.evaluate(&Fr::one()));
            let x = Fr::rand(rng);
            assert_eq!(p.evaluate(&(x * omega)), p_shifted.evaluate(&x));
        } else {
            assert_eq!(0, 1);
        }
    }
}