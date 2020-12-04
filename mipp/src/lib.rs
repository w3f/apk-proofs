use ark_bls12_381::Bls12_381;
use ark_dh_commitments::afgho16::AFGHOCommitmentG1;
use ark_dh_commitments::identity::IdentityCommitment;
use ark_dh_commitments::pedersen::PedersenCommitment;
use ark_ec::PairingEngine;
use ark_inner_products::MultiexponentiationInnerProduct;
use ark_ip_proofs::tipa::structured_scalar_message::TIPAWithSSM;
use ark_ip_proofs::tipa::TIPA;
use blake2::Blake2b;

type IP = MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>;
type GC1 = AFGHOCommitmentG1<Bls12_381>;
type SC1 = PedersenCommitment<<Bls12_381 as PairingEngine>::G1Projective>;
type IPC = IdentityCommitment<
    <Bls12_381 as PairingEngine>::G1Projective,
    <Bls12_381 as PairingEngine>::Fr,
>;

type MIPP_SRS = TIPAWithSSM<IP, GC1, IPC, Bls12_381, Blake2b>;
type TIPA_KZG = TIPA<IP, GC1, SC1, IPC, Bls12_381, Blake2b>;


#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Projective};
    use ark_dh_commitments::DoublyHomomorphicCommitment;
    use ark_dh_commitments::identity::{HomomorphicPlaceholderValue, IdentityOutput};
    use ark_ec::ProjectiveCurve;
    use ark_ff::{One, test_rng, UniformRand, Zero};
    use ark_inner_products::{InnerProduct, PairingInnerProduct};
    use rand::Rng;

    use super::*;

    #[test]
    fn mipp_srs() {
        let n = 16;

        let rng = &mut test_rng();

        let pks = (0..n)
            .map(|_| G1Projective::prime_subgroup_generator().mul(Fr::rand(rng)))
            .collect::<Vec<_>>();
        let bitmask = (0..n)
            .map(|_| if rng.gen::<bool>() { Fr::one() } else { Fr::zero() } )
            .collect::<Vec<_>>();
        let apk = IP::inner_product(&pks, &bitmask).unwrap();

        let (srs, _) = MIPP_SRS::setup(rng, n).unwrap();
        let (ck, _) = srs.get_commitment_keys();
        let v_srs = srs.get_verifier_key();

        let com_pks= PairingInnerProduct::<Bls12_381>::inner_product(&pks, &ck).unwrap();

        let proof = MIPP_SRS::prove_with_structured_scalar_message(
            &srs,
            (&pks, &bitmask),
            (&ck, &HomomorphicPlaceholderValue),
        ).unwrap();

        let proof_valid = MIPP_SRS::verify_with_unfolded_right_message(
            &v_srs,
            &HomomorphicPlaceholderValue,
            (&com_pks, &IdentityOutput(vec![apk])),
            bitmask,
            &proof,
        ).unwrap();

        assert!(proof_valid);
    }

    #[test]
    fn tipa_kzg() {
        let n = 16;

        let rng = &mut test_rng();

        let pks = (0..n)
            .map(|_| G1Projective::prime_subgroup_generator().mul(Fr::rand(rng)))
            .collect::<Vec<_>>();
        let bitmask = (0..n)
            .map(|_| if rng.gen::<bool>() { Fr::one() } else { Fr::zero() } )
            .collect::<Vec<_>>();
        let apk = IP::inner_product(&pks, &bitmask).unwrap();

        let (srs, ck_t) = TIPA_KZG::setup(rng, n).unwrap();
        let (ck_pks, ck_bitmask) = srs.get_commitment_keys();
        let v_srs = srs.get_verifier_key();

        let com_pks = GC1::commit(&ck_pks, &pks).unwrap();
        let com_bitmask = SC1::commit(&ck_bitmask, &bitmask).unwrap();
        let com_t = IPC::commit(&vec![ck_t.clone()], &vec![apk]).unwrap();

        let proof = TIPA_KZG::prove(&srs, (&pks, &bitmask), (&ck_pks, &ck_bitmask, &ck_t)).unwrap();

        let proof_valid = TIPA_KZG::verify(&v_srs, &ck_t, (&com_pks, &com_bitmask, &com_t), &proof).unwrap();

        assert!(proof_valid);
    }
}
