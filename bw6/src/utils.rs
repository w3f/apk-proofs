use ark_ff::{FftField, batch_inversion};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain, Polynomial};

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

pub fn barycentric_eval_binary_at<F: FftField>(z: F, evals: &Bitmask, domain: Radix2EvaluationDomain<F>) -> F {
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
    for b in evals.to_bits() {
        if b {
            li_inv.push(acc - F::one());
        }
        acc *= domain.group_gen_inv;
    }

    batch_inversion(&mut li_inv);
    let s: F = li_inv.iter().sum();
    z_n * s
}

pub struct LagrangeEvaluations<F: FftField> {
    pub vanishing_polynomial: F,
    pub l_0: F,
    pub l_minus_1: F,
}

pub fn lagrange_evaluations<F: FftField>(z: F, domain: Radix2EvaluationDomain<F>) -> LagrangeEvaluations<F> {
    // TODO: reuse this code with barycentric_eval methods
    let mut z_n = z; // z^n, n=2^d - domain size, so squarings only
    for _ in 0..domain.log_size_of_group {
        z_n.square_in_place();
    }

    let z_n_minus_one = z_n - F::one();
    let z_n_minus_one_div_n = z_n_minus_one * domain.size_inv;

    let mut inv = [z - F::one(), domain.group_gen * z - F::one()];
    batch_inversion(&mut inv);
    LagrangeEvaluations {
        vanishing_polynomial: z_n_minus_one,
        l_0: z_n_minus_one_div_n * inv[0],
        l_minus_1: z_n_minus_one_div_n * inv[1],
    }
}


use ark_ff::{Field, PrimeField, Zero};
use ark_ec::{AffineCurve, ProjectiveCurve};
use crate::Bitmask;

pub fn mul_then_add<G: AffineCurve>(
    bases: &[G],
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
) -> G::Projective {
    bases.iter().zip(scalars).map(|(b, s)| b.mul(*s)).sum()
}

pub fn horner<G: AffineCurve>(
    bases: &[G],
    nu: G::ScalarField,
) -> G {
    let nu = nu.into_repr();
    bases.iter().rev().fold(G::Projective::zero(), |acc, b|
        acc.mul(nu).add_mixed(b)
    ).into_affine()
}

pub fn horner_field<F: Field>(
    bases: &[F],
    nu: F,
) -> F {
    bases.iter().rev().fold(F::zero(), |acc, b| nu * acc + b)
}

/// (max_exp+1)-sized vec: 1, base, base^2,... ,base^{max_exp}
pub fn powers<F: Field>(base: F, max_exp: usize) -> Vec<F> {
    let mut result = Vec::with_capacity(max_exp + 1);
    result.push(F::one());
    if max_exp > 0 {
        result.push(base);
    }
    let mut curr = base;
    for _ in 1..max_exp {
        curr *= base;
        result.push(curr);
    };
    result
}


pub fn randomize<P, F>(
    r: F,
    polys: &[P]
) -> P
    where
        F: Field,
        P: Polynomial<F> {
    let mut res = P::zero();
    if polys.is_empty() {
        return res;
    }
    let powers = powers(r, polys.len()-1);

    powers.into_iter().zip(polys).for_each(|(r, p)| {
        res += (r, p);
    });
    res
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, One};
    use ark_poly::{Evaluations, Polynomial};
    use ark_std::{UniformRand, test_rng};
    use ark_std::convert::TryInto;
    use crate::tests::random_bits;

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

        let bitmask = Bitmask::from_bits(&random_bits(n.try_into().unwrap(), 1.0 / 2.0, rng));
        let bits_as_field_elements = bitmask.to_bits().iter()
            .map(|b| if *b { ark_bw6_761::Fr::one() } else { ark_bw6_761::Fr::zero() })
            .collect::<Vec<_>>();
        let bits_poly = Evaluations::from_vec_and_domain(bits_as_field_elements.clone(), domain).interpolate();
        let bits_poly_at_z = bits_poly.evaluate(&z);
        assert_eq!(barycentric_eval_at(z, &bits_as_field_elements, domain), bits_poly_at_z);
        assert_eq!(barycentric_eval_binary_at(z, &bitmask, domain), bits_poly_at_z);
    }

    #[test]
    pub fn test_horner() {
        let n = 10;

        let rng = &mut test_rng();

        let nu = ark_bw6_761::Fr::rand(rng);
        let bases = (0..n).map(|_| ark_bw6_761::G1Projective::rand(rng).into_affine()).collect::<Vec<_>>();

        let powers = (0..n).map(|i| nu.pow([i as u64]).into_repr()).collect::<Vec<_>>();

        assert_eq!(horner(&bases, nu), mul_then_add(&bases, &powers));
    }

    #[test]
    fn test_lagrange_evaluations() {
        let rng = &mut test_rng();

        let n = 16;
        let domain = Radix2EvaluationDomain::<ark_bw6_761::Fr>::new(n).unwrap();

        let z = ark_bw6_761::Fr::rand(rng);
        let evals = lagrange_evaluations(z, domain);
        assert_eq!(evals.vanishing_polynomial, domain.evaluate_vanishing_polynomial(z));
        let coeffs = domain.evaluate_all_lagrange_coefficients(z);
        assert_eq!(evals.l_0, coeffs[0]);
        assert_eq!(evals.l_minus_1, coeffs[n - 1]);
    }

    #[test]
    fn test_powers() {
        let rng = &mut test_rng();
        let base = ark_bw6_761::Fr::rand(rng);

        assert_eq!(powers(base, 0), vec![ark_bw6_761::Fr::one()]);
        assert_eq!(powers(base, 1), vec![ark_bw6_761::Fr::one(), base]);

        let max_exp = 255;
        let res = powers(base, max_exp);
        assert_eq!(res.len(), max_exp + 1);
        assert_eq!(res[0], ark_bw6_761::Fr::one());
        assert_eq!(res[1], base);
        assert_eq!(res[255], base.pow([255u64]));
    }
}