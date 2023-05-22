use ark_ec::CurveGroup;
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::VariableBaseMSM;
use ark_ff::{batch_inversion, Field};
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial;
use ark_std::{test_rng, UniformRand};

use crate::{final_folding_exponents, fold_points, fold_scalars, kzg};

pub struct ProverKey<E: Pairing> {
    log_m: u32,
    // m points in G1: G, tG, ..., t^(m-1)G
    ck_g1: Vec<E::G1Affine>,
    // 2m-1 points in G2: H, sH, ..., s^(2m-2)H
    ck_g2: Vec<E::G2Affine>,
    h1: E::G2Affine,
    h2: E::G2Affine,
}

pub struct VerifierKey<E: Pairing> {
    kzg_vk_g1: kzg::VerifierKeyG1<E>,
    kzg_vk_g2: kzg::VerifierKeyG2<E>,
    h1: E::G2Affine,
    h2: E::G2Affine,
}

pub struct Proof<E: Pairing> {
    // Left cross-commitments
    l_comms: Vec<PairingOutput<E>>,
    // Right cross-commitments
    r_comms: Vec<PairingOutput<E>>,
    // Final folded point
    a_final: E::G1,
    // Final folded scalar
    b_final: E::ScalarField,
    // Final commitment key in G2
    v_final: E::G2Affine,
    // Final commitment key in G1
    w_final: E::G1Affine,
    // Final commitment key argument proof for the key in G1
    kzg_proof_g1: E::G1Affine,
    // Final commitment key argument proof for the key in G2
    kzg_proof_g2: E::G2Affine,

    challenges: Challenges<E>,
}

struct Challenges<E: Pairing> {
    xs: Vec<E::ScalarField>,
    c1: E::ScalarField,
    c2: E::ScalarField,
    z: E::ScalarField,
}

impl<E: Pairing> Challenges<E> {
    fn new(log_m: usize) -> Self {
        let rng = &mut test_rng();
        Challenges {
            xs: (0..log_m).map(|_| E::ScalarField::rand(rng)).collect(),
            c1: E::ScalarField::rand(rng),
            c2: E::ScalarField::rand(rng),
            z: E::ScalarField::rand(rng),
        }
    }
}

pub fn setup<E: Pairing>(log_m: u32) -> (ProverKey<E>, VerifierKey<E>) {
    let rng = &mut test_rng();

    let m = 2usize.pow(log_m);
    let g = E::G1::rand(rng);
    let h = E::G2::rand(rng);

    let (ck_g1, kzg_vk_g1) = kzg::setup_g1::<E>(m, g, h);
    let (ck_g2, kzg_vk_g2) = kzg::setup_g2::<E>(2 * m - 1, g, h);

    let h1 = E::G2Affine::rand(rng);
    let h2 = E::G2Affine::rand(rng);

    let pk = ProverKey {
        log_m,
        ck_g1,
        ck_g2,
        h1,
        h2,
    };

    let vk = VerifierKey {
        kzg_vk_g1,
        kzg_vk_g2,
        h1,
        h2,
    };

    (pk, vk)
}

pub fn prove<E: Pairing>(pk: &ProverKey<E>, a: &[E::G1Affine], b: &[E::ScalarField]) -> Proof<E> {
    let log_m = pk.log_m as usize;
    let m = 2usize.pow(log_m as u32);

    // TODO: Fiat-Shamir
    let challenges = Challenges::<E>::new(log_m);

    let h1 = (pk.h1 * challenges.c1).into_affine();
    let h2 = (pk.h2 * challenges.c2).into_affine();

    let mut m1 = m;
    let mut a_folded = a.to_vec();
    let mut b_folded = b.to_vec();
    let mut v_folded: Vec<E::G2Affine> = pk.ck_g2.iter().cloned().step_by(2).collect();
    let mut w_folded: Vec<E::G1Affine> = pk.ck_g1.clone();
    assert_eq!(a_folded.len(), m);
    assert_eq!(b_folded.len(), m);
    assert_eq!(v_folded.len(), m);
    assert_eq!(w_folded.len(), m);

    let mut l_comms = Vec::<PairingOutput<E>>::with_capacity(log_m);
    let mut r_comms = Vec::<PairingOutput<E>>::with_capacity(log_m);

    for x in challenges.xs.iter() {
        m1 /= 2;

        let al = &a_folded[..m1];
        let ar = &a_folded[m1..];
        let bl = &b_folded[..m1];
        let br = &b_folded[m1..];
        let vl = &v_folded[..m1];
        let vr = &v_folded[m1..];
        let wl = &w_folded[..m1];
        let wr = &w_folded[m1..];

        let bl_comm = E::G1::msm(wr, bl).unwrap();
        let br_comm = E::G1::msm(wl, br).unwrap();
        let cl = E::G1::msm(ar, bl).unwrap();
        let cr = E::G1::msm(al, br).unwrap();

        // TODO: batch conversion to affine
        let bl_comm = bl_comm.into_affine();
        let br_comm = br_comm.into_affine();
        let cl = cl.into_affine();
        let cr = cr.into_affine();

        let l_vals = [ar, &[bl_comm, cl]].concat();
        let r_vals = [al, &[br_comm, cr]].concat();
        let l_leys = [vl, &[h1, h2]].concat();
        let r_keys = [vr, &[h1, h2]].concat();

        let l_comm = E::multi_pairing(l_vals, l_leys);
        let r_comm = E::multi_pairing(r_vals, r_keys);

        l_comms.push(l_comm);
        r_comms.push(r_comm);

        let x_inv = x.inverse().unwrap();

        a_folded = fold_points(al, ar, &x);
        b_folded = fold_scalars(bl, br, &x_inv);
        v_folded = fold_points(vl, vr, &x_inv);
        w_folded = fold_points(wl, wr, &x);
    }

    assert_eq!(a_folded.len(), 1);
    assert_eq!(b_folded.len(), 1);
    assert_eq!(v_folded.len(), 1);
    assert_eq!(w_folded.len(), 1);
    let a_final = a_folded[0].into();
    let b_final = b_folded[0];
    let v_final = v_folded[0];
    let w_final = w_folded[0];

    let mut xs_inv = challenges.xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let f_w = compute_final_poly_for_g1(&challenges.xs);
    let f_v = compute_final_poly_for_g2(&xs_inv);
    let kzg_proof_g1 = kzg::open_g1::<E>(&pk.ck_g1, &f_w, challenges.z);
    let kzg_proof_g2 = kzg::open_g2::<E>(&pk.ck_g2, &f_v, challenges.z);

    Proof {
        l_comms,
        r_comms,
        a_final,
        b_final,
        v_final,
        w_final,
        kzg_proof_g1,
        kzg_proof_g2,
        challenges,
    }
}

pub fn verify<E: Pairing>(vk: &VerifierKey<E>, proof: &Proof<E>, a_comm: &PairingOutput<E>, b_comm: &E::G1, c: &E::G1) {
    let challenges = &proof.challenges;

    let xs = challenges.xs.to_vec();
    let mut xs_inv = xs.clone();
    batch_inversion(xs_inv.as_mut_slice());

    let h1 = (vk.h1 * challenges.c1).into_affine();
    let h2 = (vk.h2 * challenges.c2).into_affine();

    let z = challenges.z;
    let fv_at_z = evaluate_final_poly_for_g2(&xs_inv, &z);
    let fw_at_z = evaluate_final_poly_for_g1(&xs, &z);
    assert!(kzg::verify_g2(&vk.kzg_vk_g2, proof.v_final, z, fv_at_z, proof.kzg_proof_g2));
    assert!(kzg::verify_g1(&vk.kzg_vk_g1, proof.w_final, z, fw_at_z, proof.kzg_proof_g1));

    let exps = [xs, xs_inv].concat();
    let bases = [proof.l_comms.as_slice(), proof.r_comms.as_slice()].concat();
    assert_eq!(exps.len(), bases.len());
    let comm = PairingOutput::msm(&bases, &exps).unwrap();
    // TODO: optimize pairings
    let extra = E::multi_pairing([b_comm, c], [h1, h2]);
    let comm = comm + a_comm + extra;
    let b_comm_final = proof.w_final * proof.b_final;
    let c_final = proof.a_final * proof.b_final;
    assert_eq!(comm, E::multi_pairing([proof.a_final, b_comm_final, c_final],
                                      [proof.v_final, h1, h2]));
}

// Computes the final commitment key polynomial for the contiguous SRS (in G1).
// n = 2^m, xs = [x1, ..., xm]
// fw(X) = (1 + x_m X) (1 + x_{m-1} X^2) ... (1 + x_1 X^{2^{m-1}})
// deg(fw) = n - 1
fn compute_final_poly_for_g1<F: Field>(xs: &[F]) -> DensePolynomial<F> {
    let coeffs = final_folding_exponents(xs);
    DensePolynomial::from_coefficients_vec(coeffs)
}

fn evaluate_final_poly_for_g1<F: Field>(xs: &[F], z: &F) -> F {
    evaluate_final_poly(xs, z)
}

// Computes the final commitment key polynomial for the gapped SRS (in G2).
// n = 2^m, xs = [x1, ..., xm]
// fv(X) = (1 + x_m X^2) (1 + x_{m-1} X^4) ... (1 + x_1 X^{2^m})
// deg(fv) = 2n - 2
// fv(X) = fw(X^2)
fn compute_final_poly_for_g2<F: Field>(xs: &[F]) -> DensePolynomial<F> {
    let exps = final_folding_exponents(xs);
    // interleave with 0s
    let coeffs = exps.into_iter()
        .flat_map(|exp| [exp, F::zero()])
        .collect();
    DensePolynomial::from_coefficients_vec(coeffs)
}

fn evaluate_final_poly_for_g2<F: Field>(xs: &[F], z: &F) -> F {
    let z2 = z.square();
    evaluate_final_poly(xs, &z2)
}

// Computes (1 + x_m z)(1 + x_{m-1} z^2) ... (1 + x_1 z^{2^{m-1}}).
fn evaluate_final_poly<F: Field>(xs: &[F], z: &F) -> F {
    let mut res = F::one();
    let mut z_i = z.clone();
    for x in xs.iter().rev() {
        res *= z_i * x + F::one();
        z_i = z_i.square(); //TODO: remove extra squaring
    }
    res
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
    use ark_ec::AffineRepr;
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};

    use super::*;

    fn _test_mipp<E: Pairing>() {
        let rng = &mut test_rng();

        let log_m = 2;
        let m = 2usize.pow(log_m);

        let (pk, vk) = setup::<E>(log_m);

        // Want to prove <A, b> = b1A1 + ... + bmAm = C
        let a: Vec<E::G1Affine> = (0..m).map(|_| E::G1Affine::rand(rng)).collect();
        let b: Vec<E::ScalarField> = (0..m).map(|_| E::ScalarField::rand(rng)).collect();
        let c: E::G1 = VariableBaseMSM::msm(&a, &b).unwrap();

        let v: Vec<E::G2Affine> = pk.ck_g2.iter().cloned().step_by(2).collect();
        let w: Vec<E::G1Affine> = pk.ck_g1.clone();

        // A_comm = <A, V> = e(A1, V1) * ... * e(Am, Vm)
        let a_comm: PairingOutput<E> = E::multi_pairing(&a, v);
        // b_comm = <b, W> = b1W1 * ... + bmWm
        let b_comm: E::G1 = VariableBaseMSM::msm(&w, &b).unwrap();

        let proof = prove(&pk, &a, &b);

        verify(&vk, &proof, &a_comm, &b_comm, &c);
    }

    #[test]
    fn test_mipp() {
        _test_mipp::<Bls12_381>();
    }

    fn fold_key<A: AffineRepr>(a: Vec<A>, xs: &[A::ScalarField]) -> A {
        let mut a = a;
        let mut m1 = a.len();
        for x in xs.iter() {
            m1 /= 2;
            a = fold_points(&a[..m1], &a[m1..], x);
        }
        assert_eq!(a.len(), 1);
        a[0]
    }

    #[test]
    fn test_final_ck_polynomial_for_g1() {
        let rng = &mut test_rng();

        let log_m = 8;
        let m = 2usize.pow(log_m);
        let xs: Vec<Fr> = (0..log_m).map(|_| Fr::rand(rng)).collect();

        let fw = compute_final_poly_for_g1(&xs);
        assert_eq!(fw.degree(), m - 1);

        let kzg_ck = {
            let g = G1Projective::rand(rng);
            let h = G2Projective::rand(rng);
            let (kzg_ck, _) = kzg::setup_g1::<Bls12_381>(m, g, h);
            kzg_ck
        };
        let w = kzg_ck.clone();
        let final_w = fold_key(w, &xs);
        let fw_comm = kzg::commit_g1::<Bls12_381>(&kzg_ck, &fw);
        assert_eq!(fw_comm, final_w);

        let z = Fr::rand(rng);
        let fw_at_z_1 = fw.evaluate(&z);
        let fw_at_z_2 = evaluate_final_poly_for_g1(&xs, &z);
        assert_eq!(fw_at_z_1, fw_at_z_2);
    }

    #[test]
    fn test_final_ck_polynomial_for_g2() {
        let rng = &mut test_rng();

        let log_m = 8;
        let m = 2usize.pow(log_m);
        let xs: Vec<Fr> = (0..log_m).map(|_| Fr::rand(rng)).collect();

        let fv = compute_final_poly_for_g2(&xs);
        assert_eq!(fv.degree(), 2 * m - 2);

        let kzg_ck = {
            let g = G1Projective::rand(rng);
            let h = G2Projective::rand(rng);
            let (kzg_ck, _) = kzg::setup_g2::<Bls12_381>(2 * m - 1, g, h);
            kzg_ck
        };
        let v = kzg_ck.iter().cloned().step_by(2).collect();
        let final_v = fold_key(v, &xs);
        let fv_comm = kzg::commit_g2::<Bls12_381>(&kzg_ck, &fv);
        assert_eq!(fv_comm, final_v);

        let z = Fr::rand(rng);
        let fv_at_z_1 = fv.evaluate(&z);
        let fv_at_z_2 = evaluate_final_poly_for_g2(&xs, &z);
        assert_eq!(fv_at_z_1, fv_at_z_2);
    }
}