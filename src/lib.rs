//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

#![feature(iterator_fold_self)]

pub mod bls;
pub use bls::{Signature, SecretKey, PublicKey};

use ark_ff::FftField;
use ark_std::ops::Mul;

use ark_ff::{test_rng, UniformRand, One, Zero, Field, PrimeField};
use ark_poly::{Evaluations, EvaluationDomain, GeneralEvaluationDomain, DenseOrSparsePolynomial, DensePolynomial};
use ark_poly::GeneralEvaluationDomain::Radix2;
use ark_poly_commit::kzg10::{KZG10, Powers, Randomness, VerifierKey};
use ark_poly_commit::{Polynomial, Error, PCRandomness};
use ark_bw6_761::{BW6_761, Fr};
use ark_ec::{PairingEngine, AffineCurve, ProjectiveCurve};
use ark_ec::short_weierstrass_jacobian::GroupAffine;
use ark_std::ops::{Sub};

use bitvec::vec::BitVec;
use rand::Rng;
use ark_bls12_377::{G1Projective, G1Affine};
use ark_ec::msm::VariableBaseMSM;



pub fn mul<F: Field>(s: F, p: &DensePolynomial<F>) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        p.coeffs.iter().map(|c| s * c).collect()
    )
}

pub fn prove(mut b: BitVec, pks: Vec<PublicKey>, h: G1Affine) -> Proof {
    let rng = &mut test_rng();

    let m = pks.len();

    let apk = b.iter()
        .zip(pks.iter())
        .filter(|(b, _p)| **b)
        .map(|(_b, p)| p.0)
        .sum::<G1Projective>()
        .into_affine();

    let (pks_x, pks_y): (Vec<Fr>, Vec<Fr>) = pks.iter()
        .map(|p| p.0.into_affine())
        .map(|p| (p.x, p.y))
        .unzip();



    let mut acc = vec![h;m+1];
    for (i, (b, p)) in b.iter().zip(pks.iter()).enumerate() {
        acc[i+1] = if *b {
            acc[i] + p.0.into_affine()
        } else {
            acc[i]
        }
    }

    let (mut acc_x, mut acc_y): (Vec<Fr>, Vec<Fr>) = acc.iter()
        .map(|p| (p.x, p.y))
        .unzip();

    assert_eq!(b.len(), m);
    assert_eq!(pks_x.len(), m);
    assert_eq!(pks_y.len(), m);
    assert_eq!(acc_x.len(), m+1);
    assert_eq!(acc_y.len(), m+1);
    assert_eq!(GroupAffine::new(acc_x[0], acc_y[0], false), h);
    assert_eq!(GroupAffine::new(acc_x[m], acc_y[m], false), apk + h);

    let subdomain = GeneralEvaluationDomain::<Fr>::new(m+1).unwrap();
    let n = subdomain.size();

    // Extend the computation to the whole domain
    b.resize_with(n, || false);
    // So we don't care about pks, but
    let apk_plus_h_x = acc_x[m];
    let apk_plus_h_y = acc_y[m];
    acc_x.resize_with(n, || apk_plus_h_x);
    acc_y.resize_with(n, || apk_plus_h_y);

    let b = b.iter()
        .map(|b| if *b { Fr::one() } else { Fr::zero() })
        .collect::<Vec<_>>();

    let mut acc_x_shifted = acc_x.clone();
    let mut acc_y_shifted = acc_y.clone();
    acc_x_shifted.rotate_left(1);
    acc_y_shifted.rotate_left(1);

    let b_poly = Evaluations::from_vec_and_domain(b.clone(), subdomain).interpolate();
    let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, subdomain).interpolate();
    let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, subdomain).interpolate();
    let acc_x_poly = Evaluations::from_vec_and_domain(acc_x, subdomain).interpolate();
    let acc_y_poly = Evaluations::from_vec_and_domain(acc_y, subdomain).interpolate();
    let acc_x_shifted_poly = Evaluations::from_vec_and_domain(acc_x_shifted, subdomain).interpolate();
    let acc_y_shifted_poly = Evaluations::from_vec_and_domain(acc_y_shifted, subdomain).interpolate();

    assert_eq!(b_poly.coeffs.len(), n);
    assert_eq!(b_poly.degree(), n-1);

    let domain = GeneralEvaluationDomain::<Fr>::new(4*n).unwrap();
    assert_eq!(domain.size(), 4*n);

    let B = b_poly.evaluate_over_domain_by_ref(domain);
    let x1 = acc_x_poly.evaluate_over_domain_by_ref(domain);
    let y1 = acc_y_poly.evaluate_over_domain_by_ref(domain);
    let x2 = pks_x_poly.evaluate_over_domain_by_ref(domain);
    let y2 = pks_y_poly.evaluate_over_domain_by_ref(domain);
    let x3 = acc_x_shifted_poly.evaluate_over_domain(domain);
    let y3 = acc_y_shifted_poly.evaluate_over_domain(domain);




    let (omega, omega_inv) = if let Radix2(d) = subdomain {
        (Some(d.group_gen), Some(d.group_gen_inv))
    } else {
        (None, None)
    };

    let nB = Evaluations::from_vec_and_domain(
        B.evals.iter().map(|x| Fr::one() - x).collect(),
        domain
    );


    let mut a1 =
        &(
            &B *
                &(
                    &(
                        &(
                            &(&x1 - &x2) * &(&x1 - &x2)
                        ) *
                            &(
                                &(&x1 + &x2) + &x3
                            )
                    ) -
                        &(
                            &(&y2 - &y1) * &(&y2 - &y1)
                        )
                )
        ) +
            &(
                &nB * &(&y3 - &y1)
            );

    let mut a2 =
        &(
            &B *
                &(
                    &(
                        &(&x1 - &x2) * &(&y3 + &y1)
                    ) -
                        &(
                            &(&y2 - &y1) * &(&x3 - &x1)
                        )
                )
        ) +
            &(
                &nB * &(&x3 - &x1)
            );


    let a1_poly = a1.interpolate();
    let a2_poly = a2.interpolate();
    assert_eq!(a1_poly.degree(), 4*(n-1));
    assert_eq!(a2_poly.degree(), 3*(n-1));


    let phi =  Fr::rand(rng);

    let phi_a2_poly = DensePolynomial::from_coefficients_vec(
        a2_poly.coeffs.iter().map(|c| phi * c).collect(),
    );

    let w = &a1_poly + &phi_a2_poly;

    let mut ww = Vec::with_capacity(domain.size());
    ww.push(Fr::zero());
    ww.extend(&w.coeffs);
    let mut ww = DensePolynomial::from_coefficients_vec(ww);
    ww -= &DensePolynomial::from_coefficients_vec(
        w.coeffs.iter().map(|c| omega_inv.unwrap() * c).collect(),
    );

    let (q_poly, r) = ww.divide_by_vanishing_poly(subdomain).unwrap();
    assert_eq!(r, DensePolynomial::zero());
    assert_eq!(q_poly.degree(), 3*n-3);




















    let pp = KZG10::<BW6_761>::setup(3*n-3, false, rng).unwrap();
    let (ck, vk) = KZG10::<BW6_761>::trim(&pp, 3*n-3).unwrap();
    assert_eq!(b_poly.degree(), n-1);
    assert_eq!(acc_x_poly.degree(), n-1);
    assert_eq!(acc_y_poly.degree(), n-1);
    assert_eq!(pks_x_poly.degree(), n-1);
    assert_eq!(pks_y_poly.degree(), n-1);




    let (b_comm, _) = KZG10::<BW6_761>::commit(&ck, &b_poly, None, None).unwrap();
    let (acc_x_comm, _) = KZG10::<BW6_761>::commit(&ck, &acc_x_poly, None, None).unwrap();
    let (acc_y_comm, _) = KZG10::<BW6_761>::commit(&ck, &acc_y_poly, None, None).unwrap();
    let (q_comm, _) = KZG10::<BW6_761>::commit(&ck, &q_poly, None, None).unwrap();




    let zeta  =  Fr::rand(rng);

    let proof = KZG10::<BW6_761>::open(&ck, &q_poly, zeta, &Randomness::empty()).unwrap();
    assert!(KZG10::<BW6_761>::check(&vk, &q_comm, zeta, q_poly.evaluate(zeta), &proof).is_ok());


    use crate::mul;
    let nu = ark_bw6_761::Fr::rand(rng);

    let mut powers_of_nu= vec![Fr::one()];
    let mut cur = nu;
    for _ in 0..5 {
        powers_of_nu.push(cur);
        cur *= &nu;
    }
    let constraint_polynomials1 = vec![&pks_x_poly, &pks_y_poly, &b_poly, &acc_x_poly, &acc_y_poly, &q_poly];
    assert_eq!(powers_of_nu.len(), constraint_polynomials1.len());

    let lin1_poly: DensePolynomial<Fr> = constraint_polynomials1.iter()
        .zip(powers_of_nu.clone())
        .map(|(p, c)| mul(c, p))
        .fold(DensePolynomial::zero(), |a, b| &a + &b);
    assert_eq!(lin1_poly.degree(), q_poly.degree());
    let lin1_poly_at_zeta = lin1_poly.evaluate(zeta);



    let zeta_selector = DensePolynomial::from_coefficients_vec(vec![-zeta, Fr::one()]);
    let (lin1_witness_poly, r) = DenseOrSparsePolynomial::from(&lin1_poly).divide_with_q_and_r(&zeta_selector.into()).unwrap();
    assert_eq!(r.degree(), 0);
    assert_eq!(r.coeffs[0], lin1_poly_at_zeta);

    let lin1_proof = KZG10::<BW6_761>::open(&ck, &lin1_poly, zeta, &Randomness::empty()).unwrap().w;

    let b_zeta = b_poly.evaluate(zeta);
    let pks_x_zeta = pks_x_poly.evaluate(zeta);
    let pks_y_zeta = pks_y_poly.evaluate(zeta);
    let acc_x_zeta = acc_x_poly.evaluate(zeta);
    let acc_y_zeta = acc_y_poly.evaluate(zeta);
    let acc_x_zeta_w = acc_x_poly.evaluate(zeta * omega.unwrap());
    let acc_y_zeta_w = acc_y_poly.evaluate(zeta * omega.unwrap());
    let q_zeta = q_poly.evaluate(zeta);

     Proof {
        b_comm: b_comm.0,
        acc_x_comm: acc_x_comm.0,
        acc_y_comm: acc_y_comm.0,
        q_comm: q_comm.0,

        lin1_proof,

        b_zeta,
        pks_x_zeta,
        pks_y_zeta,
        acc_x_zeta,
        acc_y_zeta,
        acc_x_zeta_w,
        acc_y_zeta_w,
        q_zeta,

        zeta,
        phi,
        nu
    }

}

pub struct Proof {
    b_comm: ark_bw6_761::G1Affine,
    acc_x_comm: ark_bw6_761::G1Affine,
    acc_y_comm: ark_bw6_761::G1Affine,
    q_comm: ark_bw6_761::G1Affine,

    lin1_proof: ark_bw6_761::G1Affine,

    pub b_zeta: Fr,
    pub pks_x_zeta: Fr,
    pub pks_y_zeta: Fr,
    pub acc_x_zeta: Fr,
    pub acc_y_zeta: Fr,
    pub acc_x_zeta_w: Fr,
    pub acc_y_zeta_w: Fr,

    pub q_zeta: Fr,

    pub zeta: ark_bw6_761::Fr,
    pub phi: ark_bw6_761::Fr,
    pub nu: ark_bw6_761::Fr,
}

pub fn verify(n: usize, vk: VerifierKey<BW6_761>, pks_x_comm: ark_bw6_761::G1Affine, pks_y_comm: ark_bw6_761::G1Affine, proof: Proof) -> bool {
    let subdomain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

    let (omega, omega_inv) = if let Radix2(d) = subdomain {
        (Some(d.group_gen), Some(d.group_gen_inv))
    } else {
        (None, None)
    };





    let nu= proof.nu;
    let mut powers_of_nu= vec![Fr::one()];
    let mut cur = nu;
    for _ in 0..5 {
        powers_of_nu.push(cur);
        cur *= &nu;
    }

    let constraint_commitments = [pks_x_comm, pks_y_comm, proof.b_comm, proof.acc_x_comm, proof.acc_y_comm, proof.q_comm];
    let constraint_evaluations1 = [proof.pks_x_zeta, proof.pks_y_zeta, proof.b_zeta, proof.acc_x_zeta, proof.acc_y_zeta, proof.q_zeta];

    let lin1_in_zeta: ark_bw6_761::Fr = constraint_evaluations1.iter()
        .zip(powers_of_nu.clone()).map(|(x, nu)| nu * x).sum();
    // assert_eq!(lin1_in_zeta, lin1_poly_at_zeta);

    let lin1_poly_committed = constraint_commitments.iter()
        .zip(powers_of_nu.clone())
        .map(|(x, nu)| x.mul(nu.into_repr()))
        .fold_first(|a, b| a + b).unwrap().into_affine();
    // assert_eq!(lin1_poly_committed, lin1_poly);


    assert!(KZG10::<BW6_761>::check(&vk, &ark_poly_commit::kzg10::Commitment(lin1_poly_committed), proof.zeta, lin1_in_zeta,
                                    &ark_poly_commit::kzg10::Proof{ w: proof.lin1_proof, random_v: None}).is_ok());


    let powers_of_nu_as_bi = powers_of_nu.iter().map(PrimeField::into_repr).collect::<Vec<_>>();

    let lin1 = VariableBaseMSM::multi_scalar_mul(
        &[proof.b_comm, proof.acc_x_comm, proof.acc_y_comm, proof.q_comm],
        &powers_of_nu_as_bi,
    );
    let lin2 = VariableBaseMSM::multi_scalar_mul(
        &[proof.acc_x_comm, proof.acc_y_comm],
        &powers_of_nu_as_bi[0..2]);

    return {
        let b = proof.b_zeta;
        let x1 = proof.acc_x_zeta;
        let y1 = proof.acc_y_zeta;
        let x2 = proof.pks_x_zeta;
        let y2 = proof.pks_y_zeta;
        let x3 = proof.acc_x_zeta_w;
        let y3 = proof.acc_y_zeta_w;

        let a1 =
            b * (
                (x1 - x2) * (x1 - x2) * (x1 + x2 + x3)
                    - (y2 - y1) * (y2 - y1)
            ) + (Fr::one() - b) * (y3 - y1);

        let a2 =
            b * (
                (x1 - x2) * (y3 + y1)
                    - (y2 - y1) * (x3 - x1)
            ) + (Fr::one() - b) * (x3 - x1);

        (a1 + proof.phi * a2) * (proof.zeta - omega_inv.unwrap()) == proof.q_zeta * subdomain.evaluate_vanishing_polynomial(proof.zeta)
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_pks_commitment() {
        let rng = &mut test_rng();
        let m = 1;
        let pks = (0..m)
            .map(|_| SecretKey::new(rng))
            .map(|sk| PublicKey::from(&sk))
            .collect::<Vec<_>>();

        let (pks_x, pks_y): (Vec<Fr>, Vec<Fr>) = pks.iter()
            .map(|p| p.0.into_affine())
            .map(|p| (p.x, p.y))
            .unzip();

        let mut b: BitVec = (0..m).map(|_| rng.gen_bool(2.0 / 3.0)).collect();
        let h = rng.gen::<G1Projective>().into_affine();

        let subdomain = GeneralEvaluationDomain::<Fr>::new(m+1).unwrap();
        let n = subdomain.size();
        let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, subdomain).interpolate();
        let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, subdomain).interpolate();


        let pp = KZG10::<BW6_761>::setup(3*n-3, false, rng).unwrap();
        let (ck, vk) = KZG10::<BW6_761>::trim(&pp, 3*n-3).unwrap();

        let (pks_x_comm, _) = KZG10::<BW6_761>::commit(&ck, &pks_x_poly, None, None).unwrap();
        let (pks_y_comm, _) = KZG10::<BW6_761>::commit(&ck, &pks_y_poly, None, None).unwrap();

        let proof = prove(b, pks, h);
        assert!(verify(n, vk,pks_x_comm.0, pks_y_comm.0, proof));
    }

    #[test]
    fn test_larger_domain() {
        let rng = &mut test_rng();
        let n = 2;
        let d1 = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let d4 = GeneralEvaluationDomain::<Fr>::new(4*n).unwrap();

        let p_evals1 = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let p_poly1 = Evaluations::from_vec_and_domain(p_evals1, d1).interpolate();

        let p_evals4 = p_poly1.evaluate_over_domain_by_ref(d4);
        let p_poly4 = p_evals4.interpolate();

        assert_eq!(p_poly1, p_poly4);
    }

    #[test]
    fn test_mul_domain() {
        let rng = &mut test_rng();
        let n = 2;
        let d1 = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let d4 = GeneralEvaluationDomain::<Fr>::new(4*n).unwrap();

        let a_evals1 = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let b_evals1 = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let a_poly1 = Evaluations::from_vec_and_domain(a_evals1, d1).interpolate();
        let b_poly1 = Evaluations::from_vec_and_domain(b_evals1, d1).interpolate();
        assert_eq!(a_poly1.degree(), n-1);
        assert_eq!(b_poly1.degree(), n-1);

        let a_evals4 = a_poly1.evaluate_over_domain_by_ref(d4);
        let b_evals4 = b_poly1.evaluate_over_domain_by_ref(d4);

        let c_evals4 = &a_evals4 * &b_evals4;
        let c_poly4 = c_evals4.interpolate();

        assert_eq!(c_poly4.degree(), 2*(n-1));
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

        if let Radix2(d) = d {
            let omega = d.group_gen;
            assert_eq!(p.evaluate(omega), p_shifted.evaluate(Fr::one()));
            let x  =  Fr::rand(rng);
            assert_eq!(p.evaluate(x * omega), p_shifted.evaluate(x));
        } else {
            assert_eq!(0, 1);
        }
    }
}