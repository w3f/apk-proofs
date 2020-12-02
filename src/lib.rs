//! Succinct proofs of a BLS public key being an aggregate key of a subset of signers given a commitment to the set of all signers' keys

// #![feature(iterator_fold_self)]

pub mod utils;
pub mod bls;
pub use bls::{Signature, SecretKey, PublicKey};

use ark_ff::{test_rng, UniformRand, One, Zero, Field, PrimeField, FftField, batch_inversion};
use ark_poly::{Evaluations, EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};
use ark_poly::Radix2EvaluationDomain;
use ark_poly_commit::kzg10::{KZG10, Powers, Randomness};
use ark_poly_commit::{PCRandomness};
use ark_ec::{AffineCurve, ProjectiveCurve, PairingEngine};
use ark_ec::short_weierstrass_jacobian::GroupAffine;


use bitvec::vec::BitVec;
use rand::Rng;
use ark_bw6_761::{BW6_761, Fr as F};
use std::time::Instant;
use ark_ec::msm::VariableBaseMSM;

type UniPoly_761 = DensePolynomial<<BW6_761 as PairingEngine>::Fr>;
type KZG_BW6 = KZG10<BW6_761, UniPoly_761>;

pub fn mul<F: Field>(s: F, p: &DensePolynomial<F>) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        p.coeffs.iter().map(|c| s * c).collect()
    )
}

pub fn mul_by_x<F: Field>(p: &DensePolynomial<F>) -> DensePolynomial<F> {
    let mut px = vec![F::zero()];
    px.extend_from_slice(&p.coeffs);
    DensePolynomial::from_coefficients_vec(px)
}

pub fn add_constant<F: FftField, D: EvaluationDomain<F>>(p: &Evaluations<F, D>, c: F, d: D) ->  Evaluations<F, D> {
    Evaluations::from_vec_and_domain(p.evals.iter().map(|x| c + x).collect(), d)
}

pub struct Params {
    domain: Radix2EvaluationDomain<F>,
    kzg_params: ark_poly_commit::kzg10::UniversalParams<BW6_761>,

    h: ark_bls12_377::G1Affine,
}

pub struct ProverKey<'a> {
    domain_size: usize,
    kzg_ck: Powers<'a, BW6_761>,
    domain: Radix2EvaluationDomain<F>,
    h: ark_bls12_377::G1Affine,
}

pub struct VerifierKey {
    domain_size: u64,
    kzg_vk: ark_poly_commit::kzg10::VerifierKey<BW6_761>,
    domain: Radix2EvaluationDomain<F>,
    h: ark_bls12_377::G1Affine,
}

pub struct LagrangeEvaluations {
    vanishing_polynomial: F,
    l_0: F,
    l_minus_1: F,
}

impl VerifierKey {
    pub fn lagrange_evaluations(&self, zeta: F) -> LagrangeEvaluations {
        let mut zeta_n = zeta;
        for _ in 0..self.domain.log_size_of_group {
            zeta_n.square_in_place();
        }
        assert_eq!(zeta_n, zeta.pow([self.domain_size]));
        let zeta_n_minus_one= zeta_n - F::one();
        let zeta_n_minus_one_div_n = zeta_n_minus_one * self.domain.size_inv;

        let mut inv = [zeta - F::one(), self.domain.group_gen * zeta - F::one()];
        batch_inversion(&mut inv);
        LagrangeEvaluations {
            vanishing_polynomial: zeta_n_minus_one,
            l_0: zeta_n_minus_one_div_n * inv[0],
            l_minus_1: zeta_n_minus_one_div_n * inv[1],
        }
    }
}

impl Params {
    pub fn new<R: Rng>(max_pks: usize, rng: &mut R) -> Self {
        let min_domain_size = max_pks + 1; // to initialize the acc with h
        let domain = Radix2EvaluationDomain::<F>::new(min_domain_size).unwrap();
        let n = domain.size();

        // deg(q) = 3n-3
        let kzg_params = KZG_BW6::setup(3*n-2, false, rng).unwrap();

        Self {
            domain,
            kzg_params,

            h: rng.gen::<ark_bls12_377::G1Projective>().into_affine(), //TODO: outside G1
        }
    }

    pub fn get_ck(&self, m: usize) -> Powers<BW6_761> {
        let powers_of_g = self.kzg_params.powers_of_g[..m].to_vec();
        Powers {
            powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
            powers_of_gamma_g: ark_std::borrow::Cow::default(),
        }
    }

    pub fn to_pk(&self) -> ProverKey {
        let n = self.domain.size();
        ProverKey {
            domain_size: n,
            kzg_ck: self.get_ck(3*n-2),
            domain: self.domain,
            h: self.h
        }
    }

    pub fn to_vk(&self) -> VerifierKey {
        let vk = ark_poly_commit::kzg10::VerifierKey {
            g: self.kzg_params.powers_of_g[0],
            gamma_g: self.kzg_params.powers_of_gamma_g[&0],
            h: self.kzg_params.h,
            beta_h: self.kzg_params.beta_h,
            prepared_h: self.kzg_params.prepared_h.clone(),
            prepared_beta_h: self.kzg_params.prepared_beta_h.clone(),
        };
        VerifierKey {
            domain_size: self.domain.size,
            kzg_vk: vk,
            domain: self.domain,
            h: self.h
        }
    }
}

struct SignerSet(Vec<PublicKey>);

impl SignerSet {
    pub fn size(&self) -> usize {
        self.0.len()
    }

    pub fn commit(&self, ck: &Powers<BW6_761>) -> (ark_bw6_761::G1Affine, ark_bw6_761::G1Affine) {
        let m = self.0.len();
        // assert_eq!(m, powers.len());
        // as now we use ifft to compute the polynomials, we require
        let domain = GeneralEvaluationDomain::<F>::new(m).unwrap();
        assert_eq!(domain.size(), ck.powers_of_g.len());

        let (pks_x, pks_y): (Vec<F>, Vec<F>) = self.0.iter()
            .map(|p| p.0.into_affine())
            .map(|p| (p.x, p.y))
            .unzip();

        let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, domain).interpolate();
        let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, domain).interpolate();

        let (pks_x_comm, _) = KZG_BW6::commit(ck, &pks_x_poly, None, None).unwrap();
        let (pks_y_comm, _) = KZG_BW6::commit(ck, &pks_y_poly, None, None).unwrap();
        (pks_x_comm.0, pks_y_comm.0)
    }

    pub fn random<R: Rng>(num_pks: usize, rng: &mut R) -> Self {
        assert!(num_pks > 1); // https://github.com/arkworks-rs/poly-commit/issues/40
        Self(
            (0..num_pks)
            .map(|_| SecretKey::new(rng))
            .map(|sk| PublicKey::from(&sk))
            .collect::<Vec<_>>()
        )
    }

    pub fn get_all(&self) -> &[PublicKey] {
        return self.0.as_slice();
    }

    pub fn get_by_mask(&self, b: &BitVec) -> Vec<&PublicKey> {
        self.0.iter().zip(b.iter()).filter(|(_p, b)| **b).map(|(p, _b)| p).collect()
    }
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

    let zeta  =  F::rand(rng);
    let zeta_omega = zeta * pk.domain.group_gen;

    let b_zeta = b_poly.evaluate(&zeta);
    let pks_x_zeta = pks_x_poly.evaluate(&zeta);
    let pks_y_zeta = pks_y_poly.evaluate(&zeta);
    let acc_x_zeta = acc_x_poly.evaluate(&zeta);
    let acc_y_zeta = acc_y_poly.evaluate(&zeta);
    let q_zeta = q_poly.evaluate(&zeta);
    let acc_x_zeta_omega = acc_x_poly.evaluate(&zeta_omega);
    let acc_y_zeta_omega = acc_y_poly.evaluate(&zeta_omega);

    let nu: u128 = rand::random();
    let nu = F::from(nu);

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

pub struct Proof {
    b_comm: ark_bw6_761::G1Affine,
    acc_x_comm: ark_bw6_761::G1Affine,
    acc_y_comm: ark_bw6_761::G1Affine,
    q_comm: ark_bw6_761::G1Affine,

    w1_proof: ark_bw6_761::G1Affine,
    w2_proof: ark_bw6_761::G1Affine,

    pub b_zeta: F,
    pub pks_x_zeta: F,
    pub pks_y_zeta: F,
    pub acc_x_zeta: F,
    pub acc_y_zeta: F,
    pub acc_x_zeta_omega: F,
    pub acc_y_zeta_omega: F,

    pub q_zeta: F,

    pub zeta: F,
    pub phi: F,
    pub nu: F,
}

pub fn verify(
    pks_x_comm: &ark_bw6_761::G1Affine,
    pks_y_comm: &ark_bw6_761::G1Affine,
    apk: PublicKey,
    bitmask: &BitVec,
    proof: &Proof,
    vk: &VerifierKey) -> bool
{
    let nu= proof.nu;

    let timer = Instant::now();
    let b_at_zeta = utils::barycentric_eval_binary_at(proof.zeta, &bitmask, vk.domain);
    println!("{}μs = accountability", timer.elapsed().as_micros());
    assert_eq!(b_at_zeta, proof.b_zeta); // accountability

    let timer = Instant::now();
    let nu_repr = nu.into_repr();
    let w2_comm = utils::horner(&[proof.acc_x_comm, proof.acc_y_comm], nu_repr).into_affine();
    let w1_comm = utils::horner(&[*pks_x_comm, *pks_y_comm, proof.b_comm, proof.q_comm, w2_comm], nu_repr).into_affine();
    println!("{}μs = multiexp", timer.elapsed().as_micros());

    let timer = Instant::now();
    let w1_zeta = utils::horner_field(&[proof.pks_x_zeta, proof.pks_y_zeta, proof.b_zeta, proof.q_zeta, proof.acc_x_zeta, proof.acc_y_zeta], nu);
    let zeta_omega = proof.zeta * vk.domain.group_gen;
    let w2_zeta_omega = utils::horner_field(&[proof.acc_x_zeta_omega, proof.acc_y_zeta_omega], nu);
    println!("{}μs = opening points evaluation", timer.elapsed().as_micros());

    let timer = Instant::now();
    let w1_comm_wrapper = ark_poly_commit::kzg10::Commitment(w1_comm);
    let w1_proof = ark_poly_commit::kzg10::Proof { w: proof.w1_proof, random_v: None };
    let w2_comm_wrapper = ark_poly_commit::kzg10::Commitment(w2_comm);
    let w2_proof = ark_poly_commit::kzg10::Proof { w: proof.w2_proof, random_v: None };
    assert!(KZG_BW6::batch_check(&vk.kzg_vk,
                                 &[w1_comm_wrapper, w2_comm_wrapper],
                                 &[proof.zeta, zeta_omega],
                                 &[w1_zeta, w2_zeta_omega],
                                 &[w1_proof, w2_proof],
                                 &mut test_rng()).unwrap());
    println!("{}μs = batched KZG openning", timer.elapsed().as_micros());

    return {
        let b = proof.b_zeta;
        let x1 = proof.acc_x_zeta;
        let y1 = proof.acc_y_zeta;
        let x2 = proof.pks_x_zeta;
        let y2 = proof.pks_y_zeta;
        let x3 = proof.acc_x_zeta_omega;
        let y3 = proof.acc_y_zeta_omega;

        let a1 =
            b * (
                (x1 - x2) * (x1 - x2) * (x1 + x2 + x3)
                    - (y2 - y1) * (y2 - y1)
            ) + (F::one() - b) * (y3 - y1);

        let a2 =
            b * (
                (x1 - x2) * (y3 + y1)
                    - (y2 - y1) * (x3 - x1)
            ) + (F::one() - b) * (x3 - x1);

        let a3 = b * (F::one() - b);

        let evals = &vk.lagrange_evaluations(proof.zeta);
        let apk = apk.0.into_affine();
        let apk_plus_h = vk.h + apk;
        let a4 = (x1 - vk.h.x) * evals.l_0 + (x1 - apk_plus_h.x) * evals.l_minus_1;
        let a5 = (y1 - vk.h.y) * evals.l_0 + (y1 - apk_plus_h.y) * evals.l_minus_1;

        let s = proof.zeta - vk.domain.group_gen_inv;
        let f = proof.phi;
        a1 * s + f * (a2 * s + f * (a3 + f * (a4 + f * a5))) == proof.q_zeta * evals.vanishing_polynomial
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use std::time::Instant;

    #[test]
    fn test_pks_commitment() {
        let num_pks = 10;

        let rng = &mut test_rng();

        let signer_set = SignerSet::random(num_pks, rng);

        let params = Params::new(signer_set.size(), rng);

        let pks_domain_size = GeneralEvaluationDomain::<F>::compute_size_of_domain(num_pks).unwrap();
        let (pks_x_comm, pks_y_comm) = signer_set.commit(&params.get_ck(pks_domain_size));

        let b: BitVec = (0..num_pks).map(|_| rng.gen_bool(2.0 / 3.0)).collect();

        let apk = bls::PublicKey::aggregate(signer_set.get_by_mask(&b));

        let proof = prove(&b, signer_set.get_all(), &params.to_pk());
        assert!(verify(&pks_x_comm, &pks_y_comm, apk, &b, &proof, &params.to_vk()));
    }

    #[test]
    fn test_bitmask() {
        use std::convert::TryInto;

        let rng = &mut test_rng();
        let n = 2u32.pow(16);
        let domain = Radix2EvaluationDomain::new(n.try_into().unwrap()).unwrap();
        let zeta = F::rand(rng);
        let b: BitVec = (0..n).map(|_| rng.gen_bool(2.0 / 3.0)).collect();
        let mut bf = b.iter()
            .map(|b| if *b { F::one() } else { F::zero() })
            .collect::<Vec<_>>();
        let b_poly = Evaluations::from_vec_and_domain(bf, domain).interpolate();
        let b_zeta = b_poly.evaluate(&zeta);

        let timer = Instant::now();
        let coeffs = domain.evaluate_all_lagrange_coefficients(zeta);
        let b_zeta2 = b.iter().zip(coeffs)
            .filter(|(b, _c)| **b)
            .map(|(_b, c)| c).sum();
        println!("{}μs - computing b(z) using domain::evaluate_all_lagrange_coefficients for n={}", timer.elapsed().as_micros(), n);

        assert_eq!(b_zeta, b_zeta2);

        let timer = Instant::now();
        let mut zeta_n = zeta;
        for _ in 0..domain.log_size_of_group {
            zeta_n.square_in_place();
        }
        zeta_n -= F::one();
        zeta_n *= &domain.size_inv;
        let mut coeffs = Vec::with_capacity(b.count_ones());
        let mut acc = zeta;
        for b in &b {
            if *b {
                coeffs.push(acc - F::one());
            }
            acc *= domain.group_gen_inv
        }
        batch_inversion(&mut coeffs);
        let sum: F = coeffs.iter().sum();
        let b_zeta3 = zeta_n * sum;
        println!("{}μs - computing b(z) using barycentric evaluation for n={}", timer.elapsed().as_micros(), n);

        assert_eq!(b_zeta, b_zeta3);
    }

    #[test]
    fn test_lagrange_evaluations() {
        let n = 16;
        let rng = &mut test_rng();
        let params = Params::new(n-1, rng);
        assert_eq!(params.domain.size(), n);

        let z = F::rand(rng);
        let evals = params.to_vk().lagrange_evaluations(z);
        assert_eq!(evals.vanishing_polynomial, params.domain.evaluate_vanishing_polynomial(z));
        let coeffs =  params.domain.evaluate_all_lagrange_coefficients(z);
        assert_eq!(evals.l_0, coeffs[0]);
        assert_eq!(evals.l_minus_1, coeffs[n-1]);
    }

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