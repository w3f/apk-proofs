use ark_ff::{PrimeField, Zero};
use ark_ec::{AffineCurve, ProjectiveCurve};

pub fn mul_then_add<G: AffineCurve>(
    bases: &[G],
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
) -> G::Projective {
    bases.iter().zip(scalars).map(|(b, s)| b.mul(*s)).sum()
}

pub fn horner<G: AffineCurve>(
    bases: &[G],
    nu: <G::ScalarField as PrimeField>::BigInt,
) -> G::Projective {
    bases.iter().rev().fold(G::Projective::zero(), |acc, b| acc.mul(nu).add_mixed(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{test_rng, UniformRand, Field};

    #[test]
    pub fn test_horner() {
        let n = 10;

        let rng = &mut test_rng();

        let nu = ark_bw6_761::Fr::rand(rng);
        let bases = (0..n).map(|_| ark_bw6_761::G1Projective::rand(rng).into_affine()).collect::<Vec<_>>();

        let powers = (0..n).map(|i| nu.pow([i as u64]).into_repr()).collect::<Vec<_>>();

        assert_eq!(horner(&bases, nu.into_repr()), mul_then_add(&bases, &powers));
    }
}