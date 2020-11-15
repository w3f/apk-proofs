use ark_ff::{FftField, batch_inversion};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};

// Evaluates a polynomial given as evaluations over a radix-2 domain (aka in Lagrange basis) in a point
// f = sum evals_i * Li
// Li(z) = (z^n -1)/n  *  w^i /(z-w^i),    li(z) := w^i /(z-w^i)
// f(z) = fi Li(z) = (z^n -1)/n sum(fi * li)
// to get li's we compute li_inv = (z-w^i)/w^i = z/w^i - 1, so we accumulate z/w^i with n multiplications,
// batch inversion costs 1 inv + 3n multiplication, and n more muls is required for li * wi,
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

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{test_rng, UniformRand};
    use ark_poly::{Evaluations, Polynomial};

    #[test]
    pub fn test_barycentric_eval_at() {
        let rng = &mut test_rng();
        let n = 2u32.pow(16);
        let domain = Radix2EvaluationDomain::new(std::convert::TryInto::try_into(n).unwrap()).unwrap();
        let evals = (0..n).map(|_| ark_bw6_761::Fr::rand(rng)).collect::<Vec<_>>();
        let z = ark_bw6_761::Fr::rand(rng);
        let poly = Evaluations::from_vec_and_domain(evals.clone(), domain).interpolate();
        let poly_at_z = poly.evaluate(&z);
        assert_eq!(barycentric_eval_at(z, &evals, domain), poly_at_z);
    }
}