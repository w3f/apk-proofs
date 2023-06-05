use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ff::Field;

pub mod ipa;
mod kzg;

fn powers<F: Field>(x: F, n: usize) -> Vec<F> {
    let x0 = F::one();
    ark_std::iter::successors(Some(x0), move |xi| Some(x * xi))
        .take(n)
        .collect()
}

fn commit_row<E: Pairing>(ck_g1: &[E::G1Affine], bit_vec: &[usize]) -> E::G1 {
    bit_vec.into_iter()
        .map(|i| ck_g1[*i])
        .sum()
}

fn commit<E: Pairing>(ck_g1: &[E::G1Affine], ck_g2: &[E::G2Affine], bit_matrix: &[Vec<usize>]) -> PairingOutput<E> {
    assert_eq!(ck_g2.len(), bit_matrix.len());
    let g1s = bit_matrix.iter().map(|bit_vec| commit_row::<E>(ck_g1, bit_vec));
    E::multi_pairing(g1s, ck_g2)
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Affine};
    use ark_ec::{CurveGroup, VariableBaseMSM};
    use ark_ff::Zero;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_std::{end_timer, start_timer, test_rng, UniformRand};
    use ark_std::rand::distributions::{Bernoulli, Standard};
    use ark_std::rand::Rng;

    use crate::ipa::{fold_points, mipp_k, mipp_u};

    use super::*;

    fn _test_matrix_vector(log_n: u32, log_m: u32) {
        let rng = &mut test_rng();

        let n = 2usize.pow(log_n);
        let m = 2usize.pow(log_m);

        let (mipp_u_pk, mipp_u_vk) = mipp_u::setup::<Bls12_381>(log_n);
        let (mipp_k_pk, mipp_k_vk) = mipp_k::setup::<Bls12_381>(log_m);

        let ck_g2_u: Vec<G2Affine> = mipp_u_pk.ck_g2.iter().cloned().step_by(2).collect();
        let ck_g2_k: Vec<G2Affine> = mipp_k_pk.ck_g2.iter().cloned().step_by(2).collect();

        let pks: Vec<G1Affine> = rng.sample_iter(Standard).take(n).collect();
        let pks_comm = Bls12_381::multi_pairing(&pks, &ck_g2_u);

        let bit_matrix: Vec<Vec<usize>> = (0..m).map(|_| {
            rng.sample_iter::<bool, _>(Bernoulli::new(2.0 / (m as f64)).unwrap())
                .take(n)
                .enumerate()
                .filter_map(|(i, b)| b.then_some(i))
                .collect()
        }).collect();

        let t_prove = start_timer!(|| format!("Proving, log(n) = {}, log(m) = {}", log_n, log_m));
        let t_bit_matrix_comm = start_timer!(|| "Committing the bit matrix");
        let row_comms: Vec<G1Projective> = bit_matrix.iter()
            .map(|bit_vec| commit_row::<Bls12_381>(&mipp_u_pk.ck_g1, bit_vec))
            .collect();
        let row_comms = G1Projective::normalize_batch(&row_comms);
        let bit_matrix_comm = Bls12_381::multi_pairing(&row_comms, &ck_g2_k);
        end_timer!(t_bit_matrix_comm);

        let t_apks_comm = start_timer!(|| "Computing and committing the apks");
        let apks: Vec<G1Projective> = bit_matrix.iter()
            .map(|bit_vec| commit_row::<Bls12_381>(&pks, bit_vec))
            .collect();
        let apks = G1Projective::normalize_batch(&apks);
        let apks_comm = Bls12_381::multi_pairing(&apks, &ck_g2_k);
        end_timer!(t_apks_comm);

        let t_mipp_k_proof = start_timer!(|| format!("MIPP-k, log(m) = {}", log_m));
        // TODO: 128-bit exps
        let lambda = Fr::rand(rng);
        let lambdas = powers(lambda, m);
        let a = G1Projective::msm(&row_comms, &lambdas).unwrap();
        let b = G1Projective::msm(&apks, &lambdas).unwrap();
        let mu = Fr::rand(rng);
        // let mipp_k_proof_a = mipp_k::prove(&mipp_k_pk, &row_comms, lambda);
        // let mipp_k_proof_b = mipp_k::prove(&mipp_k_pk, &apks, lambda);
        let comms = fold_points(&row_comms, &apks, &mu);
        let mipp_k_proof = mipp_k::prove_for_powers(&mipp_k_pk, &comms, lambda);
        end_timer!(t_mipp_k_proof);

        let t_mipp_u_proof = start_timer!(|| format!("MIPP-u, log(n) = {}", log_n));
        let mut coeffs = vec![Fr::zero(); n];
        for (row, lambda) in bit_matrix.iter().zip(lambdas) {
            for index in row {
                coeffs[*index] += lambda;
            }
        }
        let mipp_u_proof = mipp_u::prove(&mipp_u_pk, &pks, &coeffs);
        end_timer!(t_mipp_u_proof);
        end_timer!(t_prove);

        let t_verify = start_timer!(|| format!("Verification, log(n) = {}, log(m) = {}", log_n, log_m));
        // mipp_k::verify(&mipp_k_vk, &mipp_k_proof_a, &bit_matrix_comm, &lambda, &a);
        // mipp_k::verify(&mipp_k_vk, &mipp_k_proof_b, &apks_comm, &lambda, &b);
        let comm_gt = bit_matrix_comm + apks_comm * mu;
        let c = a + b * mu;
        mipp_k::verify_for_powers(&mipp_k_vk, &mipp_k_proof, &comm_gt, &lambda, &c);
        mipp_u::verify(&mipp_u_vk, &mipp_u_proof, &pks_comm, &a, &b);
        end_timer!(t_verify);

        // slot proof
        let start = 0;
        let k = 4; // commitees per slot
        let coeffs: Vec<Fr> = rng.sample_iter(Standard).take(k).collect();
        let mut r = vec![Fr::zero(); m];
        r[start..start + k].copy_from_slice(&coeffs);

        let r_row_comms = G1Projective::msm(&row_comms, &r).unwrap();
        let r_apks = G1Projective::msm(&apks, &r).unwrap();

        let mipp_k_m = mipp_k::prove_unstructured(&mipp_k_pk, &row_comms, &r);
        let mipp_k_apk = mipp_k::prove_unstructured(&mipp_k_pk, &apks, &r);
        mipp_k::verify_unstructured(&mipp_k_vk, &mipp_k_m, &bit_matrix_comm, &r, &r_row_comms);
        mipp_k::verify_unstructured(&mipp_k_vk, &mipp_k_apk, &apks_comm, &r, &r_apks);

        let mut poly = vec![Fr::zero(); n];
        for (row, r) in bit_matrix[start..start + k].iter().zip(coeffs.iter()) {
            for index in row {
                poly[*index] += r;
            }
        }
        let poly = DensePolynomial::from_coefficients_vec(poly);
        assert_eq!(kzg::commit_g1::<Bls12_381>(&mipp_u_pk.ck_g1, &poly), r_row_comms);

        let lambda = Fr::rand(rng);
        let lambdas = powers(lambda, n);

        let eval: Fr = bit_matrix[start..start + k].iter().zip(coeffs).map(|(bit_vec, r)| {
            let mut sum = Fr::zero();
            for index in bit_vec {
                sum += lambdas[*index];
            }
            sum * r
        }).sum();

        assert_eq!(poly.evaluate(&lambda), eval);
    }

    #[test]
    fn test_matrix_vector() {
        // _test_matrix_vector(10, 6);
        _test_matrix_vector(8, 5);
    }

    #[test]
    #[ignore]
    fn bench_matrix_vector() {
        _test_matrix_vector(20, 12);
    }
}