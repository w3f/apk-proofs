use ark_bw6_761;
use ark_bw6_761::Fr as F;
use ark_ec::ProjectiveCurve;
use ark_ec::short_weierstrass_jacobian::GroupAffine;
use ark_ff::{FftField, Field, One, test_rng, UniformRand, Zero};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::kzg10::{Randomness, Commitment};
use ark_poly_commit::PCRandomness;

use bitvec::vec::BitVec;

use crate::{KZG_BW6, Proof, ProverKey, PublicKey};
use ark_bw6_761::{BW6_761};

use super::ProofScheme;

//next function multiplies a polynomial (DensePolynomial<F>) by a constant s
fn mul<F: Field>(s: F, p: &DensePolynomial<F>) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        p.coeffs.iter().map(|c| s * c).collect() 
    )
}

//next function adds a constant s to a polynomial (DensePolynomial<F>)
fn add<F: Field>(s: F, p: &DensePolynomial<F>) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        p.coeffs.iter().map(|c| s * c).collect() 
    )
}

//next function multiplies a polynomial (DensePolynomial<F>) with the indeterminate X
fn mul_by_x<F: Field>(p: &DensePolynomial<F>) -> DensePolynomial<F> {
    let mut px = vec![F::zero()];
    px.extend_from_slice(&p.coeffs);
    DensePolynomial::from_coefficients_vec(px) 
}

//next function adds a constant c to all evaluations of a polynomial
fn add_constant<F: FftField, D: EvaluationDomain<F>>(p: &Evaluations<F, D>, c: F, d: D) ->  Evaluations<F, D> {
    Evaluations::from_vec_and_domain(p.evals.iter().map(|x| c + x).collect(), d)
}

//next function multiplies with a constant c all evaluations of a polynomial
fn mul_constant<F: FftField, D: EvaluationDomain<F>>(p: &Evaluations<F, D>, c: F, d: D) ->  Evaluations<F, D> {
    Evaluations::from_vec_and_domain(p.evals.iter().map(|x| c * x).collect(), d)
}

pub fn prove(b: &BitVec, pks: &[PublicKey], pk: &ProverKey, scheme: ProofScheme) -> Proof {
    let m = pks.len();

    assert_eq!(b.len(), m); //the length of the bit vector must be equal to the total number of validators
    assert!(b.count_ones() > 0); //at least one person should have signed

    let rng = &mut test_rng(); //Oana: is rng the seed for randomness generation? Alistair: is the same seed used for all randomness?

    let apk = b.iter()
        .zip(pks.iter())
        .filter(|(bit, _p)| **bit)
        .map(|(_bit, p)| p.0)
        .sum::<ark_bls12_377::G1Projective>() //TO DO
        .into_affine(); //sum of all public keys participated in signing in affine form because we are interested in field elements (x and y coordinate)

    let (pks_x, pks_y): (Vec<F>, Vec<F>) = pks.iter()
        .map(|p| p.0.into_affine())
        .map(|p| (p.x, p.y))
        .unzip(); //extracting x and y coordinate of all public keys in affine form

    let h = pk.h;
    let mut acc = vec![h;m+1];
    for (i, (b, p)) in b.iter().zip(pks.iter()).enumerate() {
        acc[i+1] = if *b {
            acc[i] + p.0.into_affine()
        } else {
            acc[i]
        }
    } //Syed writing for Alstair: what happens to the padding of acc till nearest higher 2's power? The extension happens lower there maybe it is more readable if it is done immediately.

    let (mut acc_x, mut acc_y): (Vec<F>, Vec<F>) = acc.iter()
        .map(|p| (p.x, p.y))
        .unzip(); //I think acc is a vector of affine coordinates of validator's public keys; so acc_x and acc_y are the vectors of x and y vector coordinates,respectively 

    assert_eq!(b.len(), m);
    assert_eq!(pks_x.len(), m);
    assert_eq!(pks_y.len(), m);
    assert_eq!(acc_x.len(), m+1); //the length of the vector acc_x is m+1 as we append h_x in front of acc_x.
    assert_eq!(acc_y.len(), m+1); //the length of the vector acc_y is m+1 as we append h_y in front of acc_y.
    assert_eq!(GroupAffine::new(acc_x[0], acc_y[0], false), h);
    assert_eq!(GroupAffine::new(acc_x[m], acc_y[m], false), apk + h);

    let mut b = b.iter()
        .map(|b| if *b { F::one() } else { F::zero() })
        .collect::<Vec<_>>(); //TO DO: how does the funtion collect work here?; this transforms bitvector b to 
        // corresponding vector of field elements zero and one. 

    let n = pk.domain_size; //TO DO: what exactly is pk? what is n? How big is n? Is it feasible?
    let subdomain = GeneralEvaluationDomain::<F>::new(n).unwrap();//TO DO: what is this?
    let domain = GeneralEvaluationDomain::<F>::new(4*n).unwrap(); // TO DO: what is this?
    assert_eq!(domain.size(), 4*n);
    // Extend the computation to the whole domain for polynomials higher than degree n-1

    b.resize_with(n, || F::zero());
    // So we don't care about pks, but
    let apk_plus_h_x = acc_x[m];
    let apk_plus_h_y = acc_y[m];
    acc_x.resize_with(n, || apk_plus_h_x); // resize vector acc_x to n components, filled in with component of index m
    acc_y.resize_with(n, || apk_plus_h_y);

    let mut acc_x_shifted = acc_x.clone();
    let mut acc_y_shifted = acc_y.clone();//TO DO Reminder what kind of copy this is?
    acc_x_shifted.rotate_left(1); //as per name, a rotation of vector components
    acc_y_shifted.rotate_left(1); //TODO: make sure the rotation is on correct direction

    //We need to be careful how we use l1 and ln; in my write-up and Alistair's write-up, we use l_0 and l_{n-1} and that corresponds to l_n and l_{n-1} in 
    let mut l1 = vec![F::zero(); n];
    let mut ln = vec![F::zero(); n];
    let mut lnminus1 = vec![F::zero(); n];
    l1[0] = F::one(); //create Lagrange basis vectors l1; btw, this is what I call l0 
    ln[n-1] = F::one(); //create Lagrange basis vectors ln; btw, this is what I call ln-1 


    //let b_poly = Evaluations::from_vec_and_domain(b, subdomain).interpolate();  //Initial location
    let acc_x_poly = Evaluations::from_vec_and_domain(acc_x, subdomain).interpolate();
    let acc_y_poly = Evaluations::from_vec_and_domain(acc_y, subdomain).interpolate();

    let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, subdomain).interpolate();
    let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, subdomain).interpolate();
    let acc_x_shifted_poly = Evaluations::from_vec_and_domain(acc_x_shifted, subdomain).interpolate();
    let acc_y_shifted_poly = Evaluations::from_vec_and_domain(acc_y_shifted, subdomain).interpolate();
    let l1_poly = Evaluations::from_vec_and_domain(l1, subdomain).interpolate(); //Anything related to l1 and ln are always the same and we could have done a precomputation
    let ln_poly = Evaluations::from_vec_and_domain(ln, subdomain).interpolate();
    
    //let b_comm = KZG_BW6::commit(&pk.kzg_ck, &b_poly, None, None).unwrap().0.0;  //Initial location 
    let acc_x_comm = KZG_BW6::commit(&pk.kzg_ck, &acc_x_poly, None, None).unwrap().0.0;
    let acc_y_comm = KZG_BW6::commit(&pk.kzg_ck, &acc_y_poly, None, None).unwrap().0.0;
    
    //succinct accountable variables
    let mut r_saccount;
    let mut c_saccount = vec![F::zero(); n];
    let mut c_saccount_shifted = vec![F::zero(); n];
    let mut a_saccount = vec![F::zero(); n];
    let mut a_saccount_shifted= vec![F::zero(); n];
    let mut c_saccount_poly : DensePolynomial::<F>;
    let mut a_saccount_poly : DensePolynomial::<F>;
    let mut c_saccount_shifted_poly :DensePolynomial::<F>;
    let mut a_saccount_shifted_poly :DensePolynomial::<F>;
    let mut c_saccount_comm : ark_bw6_761::G1Affine;
    let mut a_saccount_comm : ark_bw6_761::G1Affine;
    let mut acc_saccount = vec![F::zero(); n];
    let mut acc_saccount_poly: DensePolynomial::<F>;
    let mut acc_saccount_shifted = vec![F::zero(); n];
    let mut acc_saccount_shifted_poly: DensePolynomial::<F>;    
    let mut acc_saccount_comm: ark_bw6_761::G1Affine;
    let mut sum_saccount :F;
    let mut var_r_saccount :F;   
    let mut var_r_2_saccount :F;
    
    let mut a6_saccount: Evaluations<F, GeneralEvaluationDomain<F>>;// TO DO: is this the correct type? 
    let mut a7_saccount: Evaluations<F, GeneralEvaluationDomain<F>>;// TO DO: is this the correct type? 
    let mut acc_saccount_eval : Evaluations<F, GeneralEvaluationDomain<F>>;// To replace acc_saccount due to Rust move
    let mut acc_saccount_shifted_eval : Evaluations<F, GeneralEvaluationDomain<F>>;// To replace acc_saccount_shifted due to Rust move
    let mut c_saccount_eval : Evaluations<F, GeneralEvaluationDomain<F>>;// To replace c_saccount
    let mut c_saccount_shifted_eval : Evaluations<F, GeneralEvaluationDomain<F>>;// To replace c_saccount_shifted
    let mut a6_saccount_poly: DensePolynomial::<F>;
    let mut a7_saccount_poly: DensePolynomial::<F>;

    let acc_saccount_zeta: F;
    let c_saccount_zeta: F;
    let r_saccount_zeta_omega: F;
    let r_saccount_poly: DensePolynomial::<F>;

    if let ProofScheme::SuccinctAccountable = scheme {
        r_saccount = F::rand(rng); //TODO: make sure this is different than phi
	    for i in 0..n {
	        let exp1 : [u64; 1] = [(i % 256) as u64];
	        let exp2 : [u64; 1] = [(i /256) as u64];
	        c_saccount[i] = (F::one() + F::one()).pow(exp1);
            c_saccount[i] *= r_saccount.pow(exp2);
            if (i%256) == 0 {
		        a_saccount[i] = F::one();    
            } 
            else {
		        a_saccount[i] = F::zero();
            } 
            if i==0 {
                acc_saccount[i] = F::zero();
            }    
            else {
                acc_saccount[i] = acc_saccount[i-1] + b[i]*c_saccount[i];        
            }    
        }
	
        c_saccount_shifted = c_saccount.clone(); 
        c_saccount_shifted.rotate_left(1); // TO DO: make sure the rotation has correct direction

        a_saccount_shifted = a_saccount.clone(); 
        a_saccount_shifted.rotate_left(1); // TO DO: make sure the rotation has correct direction

        acc_saccount_shifted = acc_saccount.clone(); 
        acc_saccount_shifted.rotate_left(1); // TO DO: make sure the rotation has correct direction

        c_saccount_poly = Evaluations::from_vec_and_domain(c_saccount, subdomain).interpolate();
        c_saccount_shifted_poly = Evaluations::from_vec_and_domain(c_saccount_shifted, subdomain).interpolate();

        a_saccount_poly = Evaluations::from_vec_and_domain(a_saccount, subdomain).interpolate();
        a_saccount_shifted_poly = Evaluations::from_vec_and_domain(a_saccount_shifted, subdomain).interpolate();

        acc_saccount_poly = Evaluations::from_vec_and_domain(acc_saccount, subdomain).interpolate();        
        acc_saccount_shifted_poly = Evaluations::from_vec_and_domain(acc_saccount_shifted, subdomain).interpolate();

        c_saccount_comm = KZG_BW6::commit(&pk.kzg_ck, &c_saccount_poly, None, None).unwrap().0.0;
        a_saccount_comm = KZG_BW6::commit(&pk.kzg_ck, &a_saccount_poly, None, None).unwrap().0.0; //TO DO: This should be an input to prover when using succinct accountable scheme
        acc_saccount_comm = KZG_BW6::commit(&pk.kzg_ck, &acc_saccount_poly, None, None).unwrap().0.0; 

        sum_saccount = F::zero(); 

    }

    let b_poly = Evaluations::from_vec_and_domain(b, subdomain).interpolate();
    let b_comm = KZG_BW6::commit(&pk.kzg_ck, &b_poly, None, None).unwrap().0.0;
    
    assert_eq!(b_poly.coeffs.len(), n);
    assert_eq!(b_poly.degree(), n-1);// TO DO : is this always true? There are ways to get a lower degree. this assert doesn't need to hold. e.g. when b is constant

    //let domain = GeneralEvaluationDomain::<F>::new(4*n).unwrap(); //Initial location
    //assert_eq!(domain.size(), 4*n); //Initial location

    let B = b_poly.evaluate_over_domain_by_ref(domain); // TO DO: is this needed because of Rust's use of move when b? No, because of use of domain instead of subdomain
    let x1 = acc_x_poly.evaluate_over_domain_by_ref(domain);
    let y1 = acc_y_poly.evaluate_over_domain_by_ref(domain);
    let x2 = pks_x_poly.evaluate_over_domain_by_ref(domain);
    let y2 = pks_y_poly.evaluate_over_domain_by_ref(domain);
    let x3 = acc_x_shifted_poly.evaluate_over_domain(domain);
    let y3 = acc_y_shifted_poly.evaluate_over_domain(domain);// TO DO: what is the difference between evaluate_over_domain and evaluate_over_domain_by_ref ?
    let L1 = l1_poly.evaluate_over_domain(domain);
    let Ln = ln_poly.evaluate_over_domain(domain);
    //acc_saccount_shifted_eval = acc_saccount_shifted_poly.evaluate_over_domain_by_ref(domain);
    //acc_saccount_eval = acc_saccount_poly.evaluate_over_domain_by_ref(domain);
    //c_saccount_eval = c_saccount_poly.evaluate_over_domain_by_ref(domain);

    let nB = Evaluations::from_vec_and_domain(B.evals.iter().map(|x| F::one() - x).collect(),domain);

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

    let a4 = &(&acc_minus_h_x * &L1) + &(&acc_minus_h_plus_apk_x * &Ln); //a5 in Oana's writeup
    let a5 = &(&acc_minus_h_y * &L1) + &(&acc_minus_h_plus_apk_y * &Ln); //a6 in Oana's writeup
    
    //a6_saccount = &c_saccount_shifted_eval - &mul_const(&l1_poly, F::one() - var_r_saccount, domain);
    //a7_saccount = &acc_saccount_shifted_eval - &acc_saccount_eval - &c_saccount_eval * &B;

    let a1_poly = a1.interpolate();
    let a2_poly = a2.interpolate();
    let a3_poly = a3.interpolate();
    let a4_poly = a4.interpolate();
    let a5_poly = a5.interpolate();

    //a6_saccount_poly = a6_saccount.interpolate();
    //a7_saccount_poly = a7_saccount.interpolate();        

    assert_eq!(a1_poly.degree(), 4*(n-1)); // we need to assert at most not equal
    assert_eq!(a2_poly.degree(), 3*(n-1));
    assert_eq!(a3_poly.degree(), 2*(n-1));
    assert_eq!(a4_poly.degree(), 2*(n-1));
    assert_eq!(a5_poly.degree(), 2*(n-1));
    //assert_eq!(a6_poly.degree(), 2*(n-1)); // as above, the degree of a6_poly may be at most 2(n-1)
    //assert_eq!(a7_poly.degree(), 2*(n-1));

    let a1_poly_ = &mul_by_x(&a1_poly) - &mul(pk.domain.group_gen_inv, &a1_poly);
    let a2_poly_ = &mul_by_x(&a2_poly) - &mul(pk.domain.group_gen_inv, &a2_poly);

    //a6_saccount_poly_ = a6_saccount_poly - &mul(var_r_2_saccount - (F::one() + F::one()), &a_saccount_shift_poly) - &mul(1-var_r_saccount, &l1_poly);
    //a7_saccount_poly_ = a7_saccount_poly + &mul(sum, &ln_poly);
    assert_eq!(a1_poly_.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    assert_eq!(a2_poly_.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    assert_eq!(a3_poly.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    assert_eq!(a4_poly.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    assert_eq!(a5_poly.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    //assert_eq!(a6_saccount_poly_.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
    //assert_eq!(a7_saccount_poly_.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());

    let phi =  F::rand(rng); // TO DO: put a comment on what is this random is being used for? here V
    
    let mut curr = phi; //random element of the field. this is going to be replaced with Fiat Shamir
    let mut powers_of_phi = vec![curr]; //random pow
    for _ in 0..5 {
        curr *= &phi;
        powers_of_phi.push(curr);
    }
  
    let mut w = &a1_poly + &mul(powers_of_phi[0], &a2_poly); // a1 + phi a2
    w = &mul_by_x(&w) - &mul(pk.domain.group_gen_inv, &w); // X w - omega_inv w = w (X - omega_inv)
    w = &w + &mul(powers_of_phi[1], &a3_poly);
    w = &w + &mul(powers_of_phi[2], &a4_poly);
    w = &w + &mul(powers_of_phi[3], &a5_poly);
    // w = &w + &mul(powers_of_phi[4], &a6_saccount_poly);
    // w = &w + &mul(powers_of_phi[5], &a7_saccount_poly);

    let (q_poly, r) = w.divide_by_vanishing_poly(subdomain).unwrap();
    assert_eq!(r, DensePolynomial::zero());
    assert_eq!(q_poly.degree(), 3*n-3);// TO DO does this always hold?

    assert_eq!(pk.kzg_ck.powers_of_g.len(), q_poly.degree()+1);// TO DO Here we should assert only that we have enough powers in the srs
    let q_comm = KZG_BW6::commit(&pk.kzg_ck, &q_poly, None, None).unwrap().0.0;

    let zeta: F = u128::rand(rng).into(); //TO DO: why is this differently defined then phi? because we agreed on 128 bits for challenges?
    let zeta_omega = zeta * pk.domain.group_gen;

    let b_zeta = b_poly.evaluate(&zeta);
    let pks_x_zeta = pks_x_poly.evaluate(&zeta);
    let pks_y_zeta = pks_y_poly.evaluate(&zeta);
    let acc_x_zeta = acc_x_poly.evaluate(&zeta);
    let acc_y_zeta = acc_y_poly.evaluate(&zeta);
    let q_zeta = q_poly.evaluate(&zeta);//q is what I denote as t in my write-up
    let acc_x_zeta_omega = acc_x_poly.evaluate(&zeta_omega); // TO DO: these are not needed, even in the accountable and non-succinct part
    let acc_y_zeta_omega = acc_y_poly.evaluate(&zeta_omega); // TO DO: these are not needed, even in the accountable and non-succinct part
    //acc_saccount_zeta = acc_saccount_poly.evaluate(&zeta);
    //c_saccount_zeta = c_saccount_poly.evaluate(&zeta);
    //r_saccount_zeta_omega = r_saccount_poly.evaluate(&zeta_omega);
 
    let nu: F = u128::rand(rng).into(); //TO DO: why is this differently defined then phi but same as zeta?

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
