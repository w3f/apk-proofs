use ark_ip_proofs::tipa::TIPA;
use ark_inner_products::{PairingInnerProduct, ExtensionFieldElement, MultiexponentiationInnerProduct};
use ark_ip_proofs::tipa::structured_scalar_message::TIPAWithSSM;
use ark_ec::PairingEngine;
use ark_dh_commitments::afgho16::AFGHOCommitmentG1;
use ark_dh_commitments::identity::IdentityCommitment;


type MultiExpInnerProductC<P, D> = TIPAWithSSM<
    MultiexponentiationInnerProduct<<P as PairingEngine>::G1Projective>,
    AFGHOCommitmentG1<P>,
    IdentityCommitment<<P as PairingEngine>::G1Projective, <P as PairingEngine>::Fr>,
    P,
    D,
>;

#[cfg(test)]
mod tests {
    use ark_inner_products::{PairingInnerProduct, InnerProduct, MultiexponentiationInnerProduct};
    use ark_bls12_381::{Bls12_381, Fr, G1Projective};
    use ark_ip_proofs::applications::groth16_aggregation::setup_inner_product;
    use blake2::Blake2b;
    use ark_ff::{test_rng, UniformRand, Zero, One};
    use ark_ec::ProjectiveCurve;
    use bitvec::vec::BitVec;
    use ark_dh_commitments::identity::{HomomorphicPlaceholderValue, IdentityOutput};
    use rand::Rng;
    use crate::MultiExpInnerProductC;
    use ark_ip_proofs::tipa::structured_scalar_message::structured_scalar_power;

    #[test]
    fn mipp() {
        let rng = &mut test_rng();

        let n = 16;
        let sks = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let pks = sks.into_iter().map(|sk| G1Projective::prime_subgroup_generator().mul(sk)).collect::<Vec<_>>();

        let srs = setup_inner_product::<Bls12_381, Blake2b, _>(rng, n).unwrap();
        let (ck_1, ck_2) = srs.get_commitment_keys();


        let com_pk = PairingInnerProduct::<Bls12_381>::inner_product(&pks, &ck_1).unwrap();



        let bits: BitVec = (0..n).map(|_| rng.gen::<bool>()).collect();
        let bits_as_field_elements = bits.iter()
            .map(|b| if *b { Fr::one() } else { Fr::zero() })
            .collect::<Vec<_>>();

        let r = Fr::rand(rng);
        // let bits_as_field_elements = structured_scalar_power(n, &r);

        let apk = MultiexponentiationInnerProduct::<G1Projective>::inner_product(&pks, &bits_as_field_elements).unwrap();

        let tipa_proof_c = MultiExpInnerProductC::<Bls12_381, Blake2b>::prove_with_structured_scalar_message(
            &srs,
            (&pks, &bits_as_field_elements),
            (&ck_1, &HomomorphicPlaceholderValue),
        ).unwrap();

        let tipa_proof_c_valid = MultiExpInnerProductC::<Bls12_381, Blake2b>::verify_with_structured_scalar_message(
            &srs.get_verifier_key(),
            &HomomorphicPlaceholderValue,
            (&com_pk, &IdentityOutput(vec![apk])),
            &r,
            &tipa_proof_c,
        ).unwrap();

        assert!(tipa_proof_c_valid);
    }
}
