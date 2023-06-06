use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ff::{Field, Zero};
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial;
use ark_std::{end_timer, start_timer, test_rng, UniformRand};
use ark_std::rand::distributions::Standard;

use crate::ipa::{fold_points, mipp_k, mipp_u};

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

pub struct ProverKey<E: Pairing> {
    log_n: u32,
    log_m: u32,
    mipp_u_pk: mipp_u::ProverKey<E>,
    mipp_k_pk: mipp_k::ProverKey<E>,
}

pub struct VerifierKey<E: Pairing> {
    mipp_u_vk: mipp_u::VerifierKey<E>,
    mipp_k_vk: mipp_k::VerifierKey<E>,
}

pub struct EpochProof<E: Pairing> {
    mipp_u_proof: mipp_u::Proof<E>,
    mipp_k_proof: mipp_k::Proof<E>,
    a: E::G1,
    b: E::G1,
    lambda: E::ScalarField,
    mu: E::ScalarField,
}

pub struct EpochInstance<E: Pairing> {
    bit_matrix_comm: PairingOutput<E>,
    apks_comm: PairingOutput<E>,
}

pub fn setup<E: Pairing>(log_n: u32, log_m: u32) -> (ProverKey<E>, VerifierKey<E>) {
    let (mipp_u_pk, mipp_u_vk) = mipp_u::setup::<E>(log_n);
    let (mipp_k_pk, mipp_k_vk) = mipp_k::setup::<E>(log_m);
    let pk = ProverKey {
        log_n,
        log_m,
        mipp_u_pk,
        mipp_k_pk,
    };
    let vk = VerifierKey {
        mipp_u_vk,
        mipp_k_vk,
    };
    (pk, vk)
}

fn prove_epoch<E: Pairing>(pk: ProverKey<E>, params: EpochParams, pks: &[E::G1Affine], bit_matrix: Vec<Vec<usize>>) -> (EpochProof<E>, EpochInstance<E>, SlotProver<E>) {
    let rng = &mut test_rng();

    let log_n = params.log_n_signers;
    let log_m = params.log_n_slots + params.log_n_committees;

    let n = 2usize.pow(log_n);
    let m = 2usize.pow(log_m);
    let ck_g2_k: Vec<E::G2Affine> = pk.mipp_k_pk.ck_g2.iter().cloned().step_by(2).collect();

    let t_bit_matrix_comm = ark_std::start_timer!(|| "Committing the bit matrix");
    let row_comms: Vec<E::G1> = bit_matrix.iter()
        .map(|bit_row| commit_row::<E>(&pk.mipp_u_pk.ck_g1, bit_row))
        .collect();
    let row_comms = E::G1::normalize_batch(&row_comms);
    let bit_matrix_comm = E::multi_pairing(&row_comms, &ck_g2_k);
    end_timer!(t_bit_matrix_comm);

    let t_apks_comm = start_timer!(|| "Computing and committing the apks");
    let apks: Vec<E::G1> = bit_matrix.iter()
        .map(|bit_row| commit_row::<E>(&pks, bit_row))
        .collect();
    let apks = E::G1::normalize_batch(&apks);
    let apks_comm = E::multi_pairing(&apks, &ck_g2_k);
    end_timer!(t_apks_comm);

    let instance = EpochInstance {
        bit_matrix_comm,
        apks_comm,
    };

    let lambda = E::ScalarField::rand(rng);
    let lambdas = powers(lambda, m);

    let t_mipp_k_proof = start_timer!(|| format!("MIPP-k, log(m) = {}", log_m));
    let a = E::G1::msm(&row_comms, &lambdas).unwrap();
    let b = E::G1::msm(&apks, &lambdas).unwrap();
    let mu = E::ScalarField::rand(rng);
    let agg_comms = fold_points(&row_comms, &apks, &mu);
    let mipp_k_proof = mipp_k::prove_for_powers(&pk.mipp_k_pk, &agg_comms, lambda);
    end_timer!(t_mipp_k_proof);

    let mut coeffs = vec![E::ScalarField::zero(); n];
    for (bit_row, lambda) in bit_matrix.iter().zip(lambdas) {
        for index in bit_row {
            coeffs[*index] += lambda;
        }
    }

    let t_mipp_u_proof = start_timer!(|| format!("MIPP-u, log(n) = {}", log_n));
    let mipp_u_proof = mipp_u::prove(&pk.mipp_u_pk, &pks, &coeffs);
    end_timer!(t_mipp_u_proof);

    let proof = EpochProof {
        mipp_u_proof,
        mipp_k_proof,
        a,
        b,
        lambda,
        mu,
    };

    let slot_prover = SlotProver {
        epoch_params: params,
        mipp_u_pk: pk.mipp_u_pk,
        mipp_k_pk: pk.mipp_k_pk,
        bit_matrix: bit_matrix,
        row_comms,
        apks
    };

    (proof, instance, slot_prover)
}

fn verify_epoch<E: Pairing>(vk: &VerifierKey<E>,
                            pks_comm: &PairingOutput<E>,
                            instance: &EpochInstance<E>,
                            proof: &EpochProof<E>) {
    let agg_comms = instance.bit_matrix_comm + instance.apks_comm * proof.mu;
    let c = proof.a + proof.b * proof.mu;
    mipp_k::verify_for_powers(&vk.mipp_k_vk, &proof.mipp_k_proof, &agg_comms, &proof.lambda, &c);
    mipp_u::verify(&vk.mipp_u_vk, &proof.mipp_u_proof, &pks_comm, &proof.a, &proof.b);
}

#[derive(Clone)]
struct EpochParams {
    // Number of public keys.
    log_n_signers: u32,
    // Slots per epoch.
    log_n_slots: u32,
    // Committees per slot.
    log_n_committees: u32,
}

struct SlotProver<E: Pairing> {
    epoch_params: EpochParams,
    mipp_u_pk: mipp_u::ProverKey<E>,
    // TODO: mipp_k_pk is actually redundant
    mipp_k_pk: mipp_k::ProverKey<E>,
    bit_matrix: Vec<Vec<usize>>,
    row_comms: Vec<E::G1Affine>,
    apks: Vec<E::G1Affine>,
}

pub struct SlotProof<E: Pairing> {
    mipp_k_proof: mipp_k::Proof<E>,
    kzg_proof: E::G1Affine,
    r_row_comms: E::G1,
    mu: E::ScalarField,
    zeta: E::ScalarField,
}

fn prove_slot<E: Pairing>(prover: &SlotProver<E>, slot: usize, coeffs: &[E::ScalarField]) -> SlotProof<E> {
    let rng = &mut test_rng();

    let n_slots = 2usize.pow(prover.epoch_params.log_n_slots);
    let n_committees = 2usize.pow(prover.epoch_params.log_n_committees);
    let n = 2usize.pow(prover.epoch_params.log_n_signers);
    let m = n_slots * n_committees;

    assert!(slot < n_slots);
    assert_eq!(coeffs.len(), n_committees);

    let start = slot * n_committees;
    let end = start + n_committees;
    let mut r = vec![E::ScalarField::zero(); m];
    r[start..end].copy_from_slice(coeffs);

    let r_row_comms = E::G1::msm(&prover.row_comms, &r).unwrap();

    let mu = E::ScalarField::rand(rng);
    let agg_comms = fold_points(&prover.row_comms, &prover.apks, &mu);
    let mipp_k_proof = mipp_k::prove_unstructured(&prover.mipp_k_pk, &agg_comms, &r);

    let zeta = E::ScalarField::rand(rng);
    let mut poly = vec![E::ScalarField::zero(); n];
    for (bit_row, r) in prover.bit_matrix[start..end].iter().zip(coeffs) {
        for index in bit_row {
            poly[*index] += r;
        }
    }
    let poly = DensePolynomial::from_coefficients_vec(poly);
    let kzg_proof = kzg::open_g1::<E>(&prover.mipp_u_pk.ck_g1, &poly, zeta);

    SlotProof {
        mipp_k_proof,
        kzg_proof,
        r_row_comms,
        mu,
        zeta,
    }
}

pub struct SlotVerifier<E: Pairing> {
    epoch_params: EpochParams,
    mipp_u_vk: mipp_u::VerifierKey<E>,
    // TODO: mipp_k_vk is actually redundant
    mipp_k_vk: mipp_k::VerifierKey<E>,
    epoch_instance: EpochInstance<E>,

}

pub fn verify_slot<E: Pairing>(verifier: &SlotVerifier<E>, apks: &[E::G1Affine], bit_submatrix: &[Vec<usize>], coeffs: &[E::ScalarField], proof: &SlotProof<E>) {
    assert_eq!(apks.len(), coeffs.len());
    let agg_comms = verifier.epoch_instance.bit_matrix_comm + verifier.epoch_instance.apks_comm * proof.mu;
    let r_apks = E::G1::msm(apks, &coeffs).unwrap();
    let c = proof.r_row_comms + r_apks * proof.mu;
    mipp_k::verify_unstructured(&verifier.mipp_k_vk, &proof.mipp_k_proof, &agg_comms, &coeffs, &c);

    let n = 2usize.pow(verifier.epoch_params.log_n_signers);
    let zetas = powers(proof.zeta, n);
    let eval: E::ScalarField = bit_submatrix.iter().zip(coeffs).map(|(bit_row, r)| {
        let mut sum = E::ScalarField::zero();
        for index in bit_row {
            sum += zetas[*index];
        }
        sum * r
    }).sum();
    assert!(kzg::verify_g1(&verifier.mipp_u_vk.kzg_vk_g1, proof.r_row_comms.into(), proof.zeta, eval, proof.kzg_proof));
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G2Affine, G1Projective, G2Projective};
    use ark_std::{end_timer, start_timer, test_rng, UniformRand};
    use ark_std::iterable::Iterable;
    use ark_std::rand::distributions::{Bernoulli, Standard};
    use ark_std::rand::Rng;

    use super::*;

    fn _test_matrix_vector(log_n: u32, log_m: u32) {
        let rng = &mut test_rng();

        let epoch_params = EpochParams {
            log_n_signers: 8,
            log_n_slots: 2,
            log_n_committees: 3,
        };

        let log_n = epoch_params.log_n_signers;
        let log_m = epoch_params.log_n_slots + epoch_params.log_n_committees;

        let n = 2usize.pow(log_n);
        let m = 2usize.pow(log_m);
        let k = 2usize.pow(epoch_params.log_n_committees);

        let (pk, vk) = setup::<Bls12_381>(log_n, log_m);

        let ck_g2_u: Vec<G2Affine> = pk.mipp_u_pk.ck_g2.iter().cloned().step_by(2).collect();

        let g = G1Projective::rand(rng);
        let sks: Vec<Fr> = rng.sample_iter(Standard).take(n).collect();
        // let pks: Vec<G1Affine> = rng.sample_iter(Standard).take(n).collect();
        let pks: Vec<G1Projective> = sks.iter().map(|sk| g * sk).collect();
        let pks = G1Projective::normalize_batch(&pks);
        let pks_comm = Bls12_381::multi_pairing(&pks, &ck_g2_u);

        let bit_matrix: Vec<Vec<usize>> = (0..m).map(|_| {
            rng.sample_iter::<bool, _>(Bernoulli::new(2.0 / (m as f64)).unwrap())
                .take(n)
                .enumerate()
                .filter_map(|(i, b)| b.then_some(i))
                .collect()
        }).collect();

        let t_prove = start_timer!(|| format!("Proving epoch, log(n) = {}, log(m) = {}", log_n, log_m));
        let (epoch_proof, epoch_instance, slot_prover) = prove_epoch(pk, epoch_params.clone(), &pks, bit_matrix.clone());
        end_timer!(t_prove);

        let t_verify = start_timer!(|| format!("Verification, log(n) = {}, log(m) = {}", log_n, log_m));
        verify_epoch(&vk, &pks_comm, &epoch_instance, &epoch_proof);
        end_timer!(t_verify);

        // let's verify the 1st slot
        // simulate data
        let apks = &slot_prover.apks[..k];
        let messages: Vec<G2Affine> = rng.sample_iter(Standard).take(k).collect();
        let asigs: Vec<G2Projective> = bit_matrix[0..k].iter()
            .zip(messages.iter())
            .map(|(row, &m)| row.iter().map(|&index| m * sks[index]).sum())
            .collect();
        let asigs = G2Projective::normalize_batch(&asigs);
        // verify data
        // batch BLS verification equations
        let rs: Vec<Fr> = rng.sample_iter(Standard).take(k).collect();
        let rapks: Vec<G1Projective> = apks.iter()
            .zip(rs.iter())
            .map(|(&apk, r)| apk * r)
            .collect();
        let rapks = G1Projective::normalize_batch(&rapks);
        let lhs = Bls12_381::multi_pairing(rapks, messages);
        let rasigs = G2Projective::msm(&asigs, &rs).unwrap();
        let rhs = Bls12_381::pairing(g, rasigs);
        assert_eq!(lhs, rhs);


        let coeffs: Vec<Fr> = rng.sample_iter(Standard).take(k).collect();
        let slot_proof = prove_slot(&slot_prover, 0, &coeffs);

        let verifier = SlotVerifier {
            epoch_params,
            mipp_u_vk: vk.mipp_u_vk,
            mipp_k_vk: vk.mipp_k_vk,
            epoch_instance,
        };



        verify_slot(&verifier, apks, &bit_matrix[..k], &coeffs, &slot_proof);
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