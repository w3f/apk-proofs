use ark_bw6_761::Fr as F;
use ark_ec::ProjectiveCurve;
use ark_ec::short_weierstrass_jacobian::GroupAffine;
use ark_ff::{FftField, Field, One, Zero};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::kzg10::Randomness;
use ark_poly_commit::PCRandomness;
use ark_std::{UniformRand, test_rng};

use bitvec::vec::BitVec;

use crate::{KZG_BW6, Proof, ProverKey, PublicKey};

fn mul<F: Field>(s: F, p: &DensePolynomial<F>) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        p.coeffs.iter().map(|c| s * c).collect()
    )
}

fn mul_by_x<F: Field>(p: &DensePolynomial<F>) -> DensePolynomial<F> {
    let mut px = vec![F::zero()];
    px.extend_from_slice(&p.coeffs);
    DensePolynomial::from_coefficients_vec(px)
}

fn add_constant<F: FftField, D: EvaluationDomain<F>>(p: &Evaluations<F, D>, c: F, d: D) ->  Evaluations<F, D> {
    Evaluations::from_vec_and_domain(p.evals.iter().map(|x| c + x).collect(), d)
}

pub fn prove(b: &BitVec, pks: &[PublicKey], pk: &ProverKey) -> Proof {
    let m = pks.len();

    assert_eq!(b.len(), m);
    assert!(b.count_ones() > 0);

    let rng = &mut test_rng();

    let apk = b.iter()
        .zip(pks.iter())
        .filter(|(b, _p)| **b)
        .map(|(_b, p)| p.0)
        .sum::<ark_bls12_377::G1Projective>()
        .into_affine();

    let (pks_x, pks_y): (Vec<F>, Vec<F>) = pks.iter()
        .map(|p| p.0.into_affine())
        .map(|p| (p.x, p.y))
        .unzip();

    let h = pk.h;
    let mut acc = vec![h;m+1];
    for (i, (b, p)) in b.iter().zip(pks.iter()).enumerate() {
        acc[i+1] = if *b {
            acc[i] + p.0.into_affine()
        } else {
            acc[i]
        }
    }

    let (mut acc_x, mut acc_y): (Vec<F>, Vec<F>) = acc.iter()
        .map(|p| (p.x, p.y))
        .unzip();


    assert_eq!(b.len(), m);
    assert_eq!(pks_x.len(), m);
    assert_eq!(pks_y.len(), m);
    assert_eq!(acc_x.len(), m+1);
    assert_eq!(acc_y.len(), m+1);
    assert_eq!(GroupAffine::new(acc_x[0], acc_y[0], false), h);
    assert_eq!(GroupAffine::new(acc_x[m], acc_y[m], false), apk + h);

    let mut b = b.iter()
        .map(|b| if *b { F::one() } else { F::zero() })
        .collect::<Vec<_>>();

    let n = pk.domain_size;
    let subdomain = GeneralEvaluationDomain::<F>::new(n).unwrap();

    // Extend the computation to the whole domain
    b.resize_with(n, || F::zero());
    // So we don't care about pks, but
    let apk_plus_h_x = acc_x[m];
    let apk_plus_h_y = acc_y[m];
    acc_x.resize_with(n, || apk_plus_h_x);
    acc_y.resize_with(n, || apk_plus_h_y);

    let mut acc_x_shifted = acc_x.clone();
    let mut acc_y_shifted = acc_y.clone();
    acc_x_shifted.rotate_left(1);
    acc_y_shifted.rotate_left(1);

    let mut l1 = vec![F::zero(); n];
    let mut ln = vec![F::zero(); n];
    l1[0] = F::one();
    ln[n-1] = F::one();

    let b_poly = Evaluations::from_vec_and_domain(b, subdomain).interpolate();
    let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, subdomain).interpolate();
    let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, subdomain).interpolate();
    let acc_x_poly = Evaluations::from_vec_and_domain(acc_x, subdomain).interpolate();
    let acc_y_poly = Evaluations::from_vec_and_domain(acc_y, subdomain).interpolate();

    let b_comm = KZG_BW6::commit(&pk.kzg_ck, &b_poly, None, None).unwrap().0.0;
    let acc_x_comm = KZG_BW6::commit(&pk.kzg_ck, &acc_x_poly, None, None).unwrap().0.0;
    let acc_y_comm = KZG_BW6::commit(&pk.kzg_ck, &acc_y_poly, None, None).unwrap().0.0;

    let phi =  F::rand(rng);

    let acc_x_shifted_poly = Evaluations::from_vec_and_domain(acc_x_shifted, subdomain).interpolate();
    let acc_y_shifted_poly = Evaluations::from_vec_and_domain(acc_y_shifted, subdomain).interpolate();
    let l1_poly = Evaluations::from_vec_and_domain(l1, subdomain).interpolate();
    let ln_poly = Evaluations::from_vec_and_domain(ln, subdomain).interpolate();

    assert_eq!(b_poly.coeffs.len(), n);
    assert_eq!(b_poly.degree(), n-1);

    let domain = GeneralEvaluationDomain::<F>::new(4*n).unwrap();
    assert_eq!(domain.size(), 4*n);

    let B = b_poly.evaluate_over_domain_by_ref(domain);
    let x1 = acc_x_poly.evaluate_over_domain_by_ref(domain);
    let y1 = acc_y_poly.evaluate_over_domain_by_ref(domain);
    let x2 = pks_x_poly.evaluate_over_domain_by_ref(domain);
    let y2 = pks_y_poly.evaluate_over_domain_by_ref(domain);
    let x3 = acc_x_shifted_poly.evaluate_over_domain(domain);
    let y3 = acc_y_shifted_poly.evaluate_over_domain(domain);
    let L1 = l1_poly.evaluate_over_domain(domain);
    let Ln = ln_poly.evaluate_over_domain(domain);

    let nB = Evaluations::from_vec_and_domain(
        B.evals.iter().map(|x| F::one() - x).collect(),
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

    let mut a3 = &B * &nB;




    let acc_minus_h_x = add_constant(&x1, -pk.h.x, domain);
    let acc_minus_h_y = add_constant(&y1, -pk.h.y, domain);

    let acc_minus_h_plus_apk_x = add_constant(&x1, -apk_plus_h_x, domain);
    let acc_minus_h_plus_apk_y = add_constant(&y1, -apk_plus_h_y, domain);

    let a4 = &(&acc_minus_h_x * &L1) + &(&acc_minus_h_plus_apk_x * &Ln);
    let a5 = &(&acc_minus_h_y * &L1) + &(&acc_minus_h_plus_apk_y * &Ln);

    let a1_poly = a1.interpolate();
    let a2_poly = a2.interpolate();
    let a3_poly = a3.interpolate();
    let a4_poly = a4.interpolate();
    let a5_poly = a5.interpolate();

    assert_eq!(a1_poly.degree(), 4*(n-1));
    assert_eq!(a2_poly.degree(), 3*(n-1));
    assert_eq!(a3_poly.degree(), 2*(n-1));
    assert_eq!(a4_poly.degree(), 2*(n-1));
    assert_eq!(a5_poly.degree(), 2*(n-1));

    let a1_poly_ = &mul_by_x(&a1_poly) - &mul(pk.domain.group_gen_inv, &a1_poly);
    let a2_poly_ = &mul_by_x(&a2_poly) - &mul(pk.domain.group_gen_inv, &a2_poly);
    assert_eq!(a1_poly_.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    assert_eq!(a2_poly_.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    assert_eq!(a3_poly.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    assert_eq!(a4_poly.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    assert_eq!(a5_poly.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());


    let mut curr = phi;
    let mut powers_of_phi = vec![curr];
    for _ in 0..3 {
        curr *= &phi;
        powers_of_phi.push(curr);
    }

    let mut w = &a1_poly + &mul(powers_of_phi[0], &a2_poly); // a1 + phi a2
    w = &mul_by_x(&w) - &mul(pk.domain.group_gen_inv, &w); // X w - omega_inv w = w (X - omega_inv)
    w = &w + &mul(powers_of_phi[1], &a3_poly);
    w = &w + &mul(powers_of_phi[2], &a4_poly);
    w = &w + &mul(powers_of_phi[3], &a5_poly);

    let (q_poly, r) = w.divide_by_vanishing_poly(subdomain).unwrap();
    assert_eq!(r, DensePolynomial::zero());
    assert_eq!(q_poly.degree(), 3*n-3);

    assert_eq!(pk.kzg_ck.powers_of_g.len(), q_poly.degree()+1);
    let q_comm = KZG_BW6::commit(&pk.kzg_ck, &q_poly, None, None).unwrap().0.0;

    let zeta: F = u128::rand(rng).into();
    let zeta_omega = zeta * pk.domain.group_gen;

    let b_zeta = b_poly.evaluate(&zeta);
    let pks_x_zeta = pks_x_poly.evaluate(&zeta);
    let pks_y_zeta = pks_y_poly.evaluate(&zeta);
    let acc_x_zeta = acc_x_poly.evaluate(&zeta);
    let acc_y_zeta = acc_y_poly.evaluate(&zeta);
    let q_zeta = q_poly.evaluate(&zeta);
    let acc_x_zeta_omega = acc_x_poly.evaluate(&zeta_omega);
    let acc_y_zeta_omega = acc_y_poly.evaluate(&zeta_omega);

    let nu: F = u128::rand(rng).into();

    let mut curr = nu;
    let mut powers_of_nu = vec![curr];
    for _ in 0..5 {
        curr *= &nu;
        powers_of_nu.push(curr);
    }

    let w2 = &acc_x_poly + &mul(powers_of_nu[0], &acc_y_poly);
    let w2_proof = KZG_BW6::open(&pk.kzg_ck, &w2, zeta_omega, &Randomness::empty()).unwrap().w;

    let mut w1 = &pks_x_poly + &mul(powers_of_nu[0], &pks_y_poly);
    w1 = &w1 + &mul(powers_of_nu[1], &b_poly);
    w1 = &w1 + &mul(powers_of_nu[2], &q_poly);
    w1 = &w1 + &mul(powers_of_nu[3], &w2);
    let w1_proof = KZG_BW6::open(&pk.kzg_ck, &w1, zeta, &Randomness::empty()).unwrap().w;

    Proof {
        b_comm,
        acc_x_comm,
        acc_y_comm,
        q_comm,

        w1_proof,
        w2_proof,

        b_zeta,
        pks_x_zeta,
        pks_y_zeta,
        acc_x_zeta,
        acc_y_zeta,
        acc_x_zeta_omega,
        acc_y_zeta_omega,
        q_zeta,

        zeta,
        phi,
        nu,
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;

    #[test]
    fn test_larger_domain() {
        let rng = &mut test_rng();
        let n = 2;
        let d1 = GeneralEvaluationDomain::<F>::new(n).unwrap();
        let d4 = GeneralEvaluationDomain::<F>::new(4*n).unwrap();

        let p_evals1 = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let p_poly1 = Evaluations::from_vec_and_domain(p_evals1, d1).interpolate();

        let p_evals4 = p_poly1.evaluate_over_domain_by_ref(d4);
        let p_poly4 = p_evals4.interpolate();

        assert_eq!(p_poly1, p_poly4);
    }

    #[test]
    fn test_mul_domain() {
        let rng = &mut test_rng();
        let n = 2;
        let d1 = GeneralEvaluationDomain::<F>::new(n).unwrap();
        let d4 = GeneralEvaluationDomain::<F>::new(4*n).unwrap();

        let a_evals1 = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let b_evals1 = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
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
        let d = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let p_evals = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let mut p_evals_shifted = p_evals.clone();
        p_evals_shifted.rotate_left(1);

        let p = Evaluations::from_vec_and_domain(p_evals, d).interpolate();
        let p_shifted = Evaluations::from_vec_and_domain(p_evals_shifted, d).interpolate();

        if let ark_poly::GeneralEvaluationDomain::Radix2(d) = d {
            let omega = d.group_gen;
            assert_eq!(p.evaluate(&omega), p_shifted.evaluate(&F::one()));
            let x  =  F::rand(rng);
            assert_eq!(p.evaluate(&(x * omega)), p_shifted.evaluate(&x));
        } else {
            assert_eq!(0, 1);
        }
    }
}