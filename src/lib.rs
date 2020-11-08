//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

pub mod bls;
pub use bls::{Signature, SecretKey, PublicKey};

#[cfg(test)]
pub mod tests {
    use ark_ff::{test_rng, UniformRand, One, Zero, Field};
    use ark_poly::{Evaluations, EvaluationDomain, GeneralEvaluationDomain, DensePolynomial};
    use ark_poly::GeneralEvaluationDomain::Radix2;
    use ark_poly_commit::kzg10::{KZG10, Powers};
    use ark_poly_commit::{Polynomial, Error};
    use ark_bw6_761::{BW6_761, Fr};
    use ark_ec::{PairingEngine, ProjectiveCurve};
    use ark_ec::short_weierstrass_jacobian::GroupAffine;

    use ark_std::ops::{Sub};



    use crate::{SecretKey, PublicKey};

    use bitvec::vec::BitVec;
    use rand::Rng;
    use ark_bls12_377::{G1Projective, G1Affine};

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

        println!("{}/{} bits set", b.count_ones(), b.len());

        let h= rng.gen::<G1Projective>().into_affine();

        let mut acc = vec![h;m+1];
        for (i, (b, p)) in b.iter().zip(pks.iter()).enumerate() {
            acc[i+1] = if *b {
                acc[i] + p.0.into_affine()
            } else {
                acc[i]
            }
        }

        let apk = b.iter()
            .zip(pks.iter())
            .filter(|(b, _p)| **b)
            .map(|(_b, p)| p.0)
            .sum::<G1Projective>()
            .into_affine();

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

        let (q, r) = ww.divide_by_vanishing_poly(subdomain).unwrap();
        assert_eq!(r, DensePolynomial::zero());
        assert_eq!(q.degree(), 3*n-3);

        // let pp = KZG10::<BW6_761>::setup(n, false, rng).unwrap();
        //
        // let powers = Powers {
        //     powers_of_g: ark_std::borrow::Cow::Owned(pp.powers_of_g),
        //     powers_of_gamma_g: ark_std::borrow::Cow::default(),
        // };

        // let (comm, _) = KZG10::<BW6_761>::commit(&powers, &p, None, None).unwrap();

        let zeta  =  Fr::rand(rng);

        struct Proof {
            pub b_zeta: Fr,
            pub pks_x_zeta: Fr,
            pub pks_y_zeta: Fr,
            pub acc_x_zeta: Fr,
            pub acc_y_zeta: Fr,
            pub acc_x_zeta_w: Fr,
            pub acc_y_zeta_w: Fr,

            pub q_zeta: Fr
        }

        let b_zeta = b_poly.evaluate(zeta);
        let pks_x_zeta = pks_x_poly.evaluate(zeta);
        let pks_y_zeta = pks_y_poly.evaluate(zeta);
        let acc_x_zeta = acc_x_poly.evaluate(zeta);
        let acc_y_zeta = acc_y_poly.evaluate(zeta);
        let acc_x_zeta_w = acc_x_poly.evaluate(zeta * omega.unwrap());
        let acc_y_zeta_w = acc_y_poly.evaluate(zeta * omega.unwrap());
        let q_zeta = q.evaluate(zeta);

        let proof = Proof {
            b_zeta,
            pks_x_zeta,
            pks_y_zeta,
            acc_x_zeta,
            acc_y_zeta,
            acc_x_zeta_w,
            acc_y_zeta_w,
            q_zeta
        };

        // Verifier needs to evaluate a1
        {
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

            assert_eq!((a1 + phi * a2) * (zeta - omega_inv.unwrap()), proof.q_zeta * subdomain.evaluate_vanishing_polynomial(zeta));
        }
    }
}