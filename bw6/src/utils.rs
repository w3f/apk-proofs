use ark_ff::{FftField, batch_inversion};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};

// Evaluates a polynomial represented as evaluations over a radix-2 domain (aka in Lagrange basis) at a point.
// f = sum(fi * Li), where Li is the i-th Lagrange basis polynomial, and fi = f(w^i)
// Li(z) = (z^n -1)/n * li(z), where li(z) := w^i /(z-w^i), see https://hackmd.io/xTta-c--SFyOv9Kl3Q9Jjw
// Then f(z) = sum(fi * Li(z)) = (z^n -1)/n sum(fi * li(z))
// To avoid inversions when computing li(z) we instead compute li_inv = (z-w^i)/w^i = z/w^i - 1 = z * w_inv^i - 1
// using n multiplications to accumulate z * w_inv^i and then perform batch inversion.
// Batch inversion costs 1 inv + 3n muls, and n more muls is required for fi * li(z),
// resulting in around 1 inv + 5n muls
pub fn barycentric_eval_at<F: FftField>(z: F, evals: &Vec<F>, domain: Radix2EvaluationDomain<F>) -> F {
    let n = domain.size();
    assert_eq!(evals.len(), n);
    // let timer_z_n =  std::time::Instant::now();
    let mut z_n = z; // z^n, n=2^d - domain size, so squarings only
    for _ in 0..domain.log_size_of_group {
        z_n.square_in_place();
    }
    // println!("{}μs z^n for log_n={}", timer_z_n.elapsed().as_micros(), domain.log_size_of_group);
    z_n -= F::one();
    z_n *= &domain.size_inv; // (z^n-1)/n

    let mut li_inv = Vec::with_capacity(n);
    let mut acc = z;
    for _ in 0..n {
        li_inv.push(acc - F::one());
        acc *= domain.group_gen_inv;
    }
    batch_inversion(&mut li_inv);
    let s = evals.iter().zip(li_inv).map(|(fi, li)| li * fi).sum::<F>();
    z_n * s
}

pub fn barycentric_eval_binary_at<F: FftField>(z: F, evals: &BitVec, domain: Radix2EvaluationDomain<F>) -> F {
    // let timer_z_n =  std::time::Instant::now();
    let mut z_n = z; // z^n, n=2^d - domain size, so squarings only
    for _ in 0..domain.log_size_of_group {
        z_n.square_in_place();
    }
    // println!("{}μs z^n for log_n={}", timer_z_n.elapsed().as_micros(), domain.log_size_of_group);
    z_n -= F::one();
    z_n *= &domain.size_inv; // (z^n-1)/n

    let mut li_inv = Vec::with_capacity(evals.count_ones());
    let mut acc = z;
    for b in evals {
        if *b {
            li_inv.push(acc - F::one());
        }
        acc *= domain.group_gen_inv;
    }

    batch_inversion(&mut li_inv);
    let s: F = li_inv.iter().sum();
    z_n * s
}


use ark_ff::{Field, PrimeField, Zero};
use ark_ec::{AffineCurve, ProjectiveCurve};
use bitvec::vec::BitVec;

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

pub fn horner_field<F: Field>(
    bases: &[F],
    nu: F,
) -> F {
    bases.iter().rev().fold(F::zero(), |acc, b| nu * acc + b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{test_rng, UniformRand, Field, One};
    use ark_poly::{Evaluations, Polynomial};
    use rand::Rng;

    #[test]
    pub fn test_barycentric_eval() {
        let rng = &mut test_rng();
        let n = 2u32.pow(16);
        let domain = Radix2EvaluationDomain::new(std::convert::TryInto::try_into(n).unwrap()).unwrap();
        let z = ark_bw6_761::Fr::rand(rng);

        let evals = (0..n).map(|_| ark_bw6_761::Fr::rand(rng)).collect::<Vec<_>>();
        let poly = Evaluations::from_vec_and_domain(evals.clone(), domain).interpolate();
        let poly_at_z = poly.evaluate(&z);
        assert_eq!(barycentric_eval_at(z, &evals, domain), poly_at_z);

        let bits: BitVec = (0..n).map(|_| rng.gen::<bool>()).collect();
        let bits_as_field_elements = bits.iter()
            .map(|b| if *b { ark_bw6_761::Fr::one() } else { ark_bw6_761::Fr::zero() })
            .collect::<Vec<_>>();
        let bits_poly = Evaluations::from_vec_and_domain(bits_as_field_elements.clone(), domain).interpolate();
        let bits_poly_at_z = bits_poly.evaluate(&z);
        assert_eq!(barycentric_eval_at(z, &bits_as_field_elements, domain), bits_poly_at_z);
        assert_eq!(barycentric_eval_binary_at(z, &bits, domain), bits_poly_at_z);
    }

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