use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::VariableBaseMSM;
use ark_ff::{batch_inversion, Field};
use ark_std::{end_timer, start_timer, test_rng};
use ark_std::rand::Rng;

use crate::ipa::{final_folding_exponents, fold_points};

pub struct Proof<E: Pairing> {
    // Left cross-commitments
    l_comms: Vec<PairingOutput<E>>,
    // Right cross-commitments
    r_comms: Vec<PairingOutput<E>>,
    // Challenges
    xs: Vec<E::ScalarField>,
}

pub fn prove<E: Pairing>(log_n: usize, a: &[E::G1Affine], b: &[E::G2Affine]) -> Proof<E> {
    let rng = &mut test_rng();

    let n = 2usize.pow(log_n as u32);

    let xs: Vec<E::ScalarField> = (0..log_n).map(|_| E::ScalarField::from(rng.gen::<u128>())).collect();

    let mut n1 = n;
    let mut a_folded = a.to_vec();
    let mut b_folded = b.to_vec();
    assert_eq!(a_folded.len(), n);
    assert_eq!(b_folded.len(), n);

    let mut l_comms = Vec::<PairingOutput<E>>::with_capacity(log_n);
    let mut r_comms = Vec::<PairingOutput<E>>::with_capacity(log_n);

    let t_gipa = start_timer!(|| "GIPA");
    for x in xs.iter() {
        let t_round = start_timer!(|| format!("ROUND: n = {}", n1));
        n1 /= 2;

        let al = &a_folded[..n1];
        let ar = &a_folded[n1..];
        let bl = &b_folded[..n1];
        let br = &b_folded[n1..];

        let t_multipairing = start_timer!(|| format!("2 x {}-multipairing", n1));
        let l_comm = E::multi_pairing(ar, bl);
        let r_comm = E::multi_pairing(al, br);
        end_timer!(t_multipairing);
        l_comms.push(l_comm);
        r_comms.push(r_comm);

        let x_inv = x.inverse().unwrap();
        let t_folding = start_timer!(|| format!("{}-folding in G1 and G2", n1));
        a_folded = fold_points(al, ar, &x);
        b_folded = fold_points(bl, br, &x_inv);
        end_timer!(t_folding);
        end_timer!(t_round);
    }
    end_timer!(t_gipa);

    assert_eq!(a_folded.len(), 1);
    assert_eq!(b_folded.len(), 1);

    Proof {
        l_comms,
        r_comms,
        xs,
    }
}

pub fn verify<E: Pairing>(proof: &Proof<E>, a: &[E::G1Affine], b: &[E::G2Affine], c: &PairingOutput<E>) {
    let t_inv = start_timer!(|| "Challenges' inverses");
    let mut xs_inv = proof.xs.clone();
    batch_inversion(xs_inv.as_mut_slice());
    end_timer!(t_inv);

    let bases = [proof.l_comms.as_slice(), proof.r_comms.as_slice()].concat();
    let scalars = [proof.xs.as_slice(), xs_inv.as_slice()].concat();
    let c = c + PairingOutput::<E>::msm(&bases, &scalars).unwrap();

    let t_exps = start_timer!(|| "Final exponents");
    let exps = final_folding_exponents(&proof.xs);
    let exps_inv = final_folding_exponents(&xs_inv);
    end_timer!(t_exps);

    let t_msm = start_timer!(|| "MSMs");
    let a1 = E::G1::msm(&a, &exps).unwrap();
    let b1 = E::G2::msm(&b, &exps_inv).unwrap();
    end_timer!(t_msm);

    let t_pairing = start_timer!(|| "Final pairing");
    let c1 = E::pairing(&a1, &b1);
    end_timer!(t_pairing);

    assert_eq!(c, c1);
}


#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_std::{test_rng, UniformRand};

    use super::*;

    fn _test_sipp<E: Pairing>() {
        let rng = &mut test_rng();

        let log_n = 8;
        let n = 2usize.pow(log_n as u32);

        // Want to prove <A, B> = e(A1, B1) * ... * e(An, Bn) = C
        let a: Vec<E::G1Affine> = (0..n).map(|_| E::G1Affine::rand(rng)).collect();
        let b: Vec<E::G2Affine> = (0..n).map(|_| E::G2Affine::rand(rng)).collect();
        let c = E::multi_pairing(&a, &b);

        let t_prove = start_timer!(|| format!("SIPP, log(n) = {}", log_n));
        let proof = prove(log_n, &a, &b);
        end_timer!(t_prove);

        verify(&proof, &a, &b, &c);
    }

    #[test]
    fn test_sipp() {
        _test_sipp::<Bls12_381>();
    }
}