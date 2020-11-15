use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ark_ff::{Field, test_rng, UniformRand};

extern crate apk_proofs;

fn barycentric_evaluation<F: Field>(c: &mut Criterion, n: u32) {
    use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};

    let rng = &mut test_rng();
    let n = std::convert::TryInto::try_into(n).unwrap();
    let domain = Radix2EvaluationDomain::new(n).unwrap();
    let evals = (0..n).map(|_| ark_bw6_761::Fr::rand(rng)).collect::<Vec<_>>();
    let z = ark_bw6_761::Fr::rand(rng);
    c.bench_function("barycentric_evaluation", move |b| {
        b.iter(|| {
            apk_proofs::utils::barycentric_eval_at(z, &evals, domain)
        })
    });
}

fn criterion_benchmark(c: &mut Criterion) {
    barycentric_evaluation::<ark_bw6_761::Fr>(c, 2u32.pow(16));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);