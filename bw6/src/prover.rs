use ark_bw6_761::{BW6_761, Fr};
use ark_ec::ProjectiveCurve;
use ark_ff::{One, Zero};
use ark_poly::{Evaluations, Polynomial, Radix2EvaluationDomain};
use ark_poly::univariate::DensePolynomial;
use merlin::Transcript;

use crate::{KZG_BW6, Proof, point_in_g1_complement, Bitmask, BasicRegisterCommitments, RegisterCommitments, ExtendedRegisterCommitments, PackedRegisterCommitments, AdditionalCommitments};
use crate::transcript::ApkTranscript;
use crate::signer_set::SignerSetCommitment;
use crate::kzg::ProverKey;
use crate::bls::PublicKey;
use crate::domains::Domains;
use crate::constraints::{Registers, RegisterEvaluations, SuccinctAccountableRegisterEvaluations, SuccinctlyAccountableRegisters, BasicRegisterEvaluations};
use crate::piop::{PiopDecorator, AdditionalRegisterPolynomials, PackedAccountabilityRegisterPolynomials};


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

    pub fn prove_simple(&self, bitmask: &Bitmask) -> Proof<BasicRegisterEvaluations, BasicRegisterCommitments> {
        self.prove::<
            BasicRegisterEvaluations,
            (),
            (),
            BasicRegisterCommitments,
            Registers
        >(bitmask)
    }

    pub fn prove_packed(&self, bitmask: &Bitmask) -> Proof<SuccinctAccountableRegisterEvaluations, ExtendedRegisterCommitments<PackedRegisterCommitments>> {
        self.prove::<
            SuccinctAccountableRegisterEvaluations,
            PackedRegisterCommitments,
            PackedAccountabilityRegisterPolynomials,
            ExtendedRegisterCommitments<PackedRegisterCommitments>,
            SuccinctlyAccountableRegisters
        >(bitmask)
    }

    #[allow(non_snake_case)]
    fn prove<E, AC, AP, C, D>(&self, bitmask: &Bitmask) -> Proof<E, C>
    where
        E: RegisterEvaluations,
        AC: AdditionalCommitments,
        AP: AdditionalRegisterPolynomials<AC>,
        C: RegisterCommitments<AC>,
        D: PiopDecorator<E, AP>,
    {
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

        // 1. Compute and commit to the basic registers.
        let registers = Registers::new(self.domains.clone(), bitmask, pks);
        let b_poly = registers.get_bitmask_register_polynomial();
        let (acc_x_poly, acc_y_poly) = registers.get_partial_sums_register_polynomials();
        let b_comm = KZG_BW6::commit(&self.params.kzg_pk, &b_poly);
        let acc_x_comm = KZG_BW6::commit(&self.params.kzg_pk, &acc_x_poly);
        let acc_y_comm = KZG_BW6::commit(&self.params.kzg_pk, &acc_y_poly);
        let basic_commitments = BasicRegisterCommitments {
            b_comm,
            acc_comm: (acc_x_comm, acc_y_comm)
        };
        transcript.append_basic_commitments(&basic_commitments);

        // 2. Receive bitmask aggregation challenge,
        // compute and commit to succinct accountability registers.
        let r = transcript.get_128_bit_challenge(b"r"); // bitmask aggregation challenge
        let acc_registers = D::wrap(registers, b, r);
        let acc_register_polynomials = acc_registers.get_accountable_register_polynomials();
        let acc_register_commitments = acc_register_polynomials.map(|polys| {
            polys.commit(|p| KZG_BW6::commit(&self.params.kzg_pk, &p))
        });
        let register_commitments = C::wrap(basic_commitments, acc_register_commitments);
        transcript.append_accountability_commitments(register_commitments.get_additional_commitments());

        // 3. Receive constraint aggregation challenge,
        // compute and commit to the quotient polynomial.
        let phi = transcript.get_128_bit_challenge(b"phi"); // constraint polynomials batching challenge
        let q_poly = acc_registers.compute_quotient_polynomial(phi, &self.domains);
        assert_eq!(self.params.kzg_pk.max_degree(), q_poly.degree()); //TODO: check at the prover creation
        assert_eq!(q_poly.degree(), 3 * n - 3);
        let q_comm = KZG_BW6::commit(&self.params.kzg_pk, &q_poly);
        transcript.append_proof_point(b"q_comm", &q_comm);

        // 4. Receive the evaluation point,
        // evaluate register polynomials and the quotient polynomial,
        // and commit to the evaluations.
        let zeta = transcript.get_128_bit_challenge(b"zeta"); // evaluation point challenge
        let register_evaluations = acc_registers.evaluate_register_polynomials(zeta);
        let q_zeta = q_poly.evaluate(&zeta);
        transcript.append_evals(&register_evaluations);
        transcript.append_proof_scalar(b"q_zeta", &q_zeta);


        // 5. Compute the linearization polynomial,
        // evaluate it at the shifted evaluation point,
        // and commit to the evaluation.
        let zeta_omega = zeta * self.domains.omega;
        let zeta_minus_omega_inv = zeta - self.domains.omega_inv;
        let r_poly = acc_registers.compute_linearization_polynomial(&register_evaluations, phi, zeta_minus_omega_inv);
        let r_zeta_omega = r_poly.evaluate(&zeta_omega);
        transcript.append_proof_scalar(b"r_zeta_omega", &r_zeta_omega);

        // 6. Receive the polynomials aggregation challenge,
        // open the aggregated polynomial at the evaluation point,
        // and the linearization polynomial at the shifted evaluation point,
        // and commit to the opening proofs.
        let nu: Fr = transcript.get_128_bit_challenge(b"nu"); // KZG opening batching challenge
        let mut register_polynomials = acc_registers.get_all_register_polynomials();
        register_polynomials.push(q_poly);
        let w_poly = KZG_BW6::aggregate_polynomials(nu, &register_polynomials);
        let w_at_zeta_proof = KZG_BW6::open(&self.params.kzg_pk, &w_poly, zeta);
        transcript.append_proof_point(b"w_at_zeta_proof", &w_at_zeta_proof);

        let r_at_zeta_omega_proof = KZG_BW6::open(&self.params.kzg_pk, &r_poly, zeta_omega);
        transcript.append_proof_point(b"r_at_zeta_omega_proof", &r_at_zeta_omega_proof);

        // Finally, compose the proof.
        Proof {
            register_commitments,
            // phi <-
            q_comm,
            // zeta <-
            register_evaluations,
            q_zeta,
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