use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ark_ff::{Field, PrimeField, Zero, test_rng, UniformRand};
use ark_ec::{AffineCurve, ProjectiveCurve};

extern crate apk_proofs;


fn mul_then_add<G: AffineCurve>(
    bases: &[G],
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
) -> G::Projective {
    bases.iter().zip(scalars).map(|(b, s)| b.mul(*s)).sum()
}

fn msm<G: AffineCurve>(c: &mut Criterion, n: usize) {

    use ark_ec::msm::VariableBaseMSM;

    let rng = &mut test_rng();

    let scalars = (0..n)
        .map(|_| G::ScalarField::rand(rng).into_repr())
        .collect::<Vec<_>>();
    let bases = (0..n)
        .map(|_| G::Projective::rand(rng).into_affine())
        .collect::<Vec<_>>();

    let scalars_clone = scalars.clone();
    let bases_clone = bases.clone();
    c.bench_function("ark_ec::msm::VariableBaseMSM", move |b| {
        b.iter(|| VariableBaseMSM::multi_scalar_mul(black_box(&bases_clone), black_box(&scalars_clone)))
    });

    let scalars_clone = scalars.clone();
    let bases_clone = bases.clone();
    c.bench_function("naive mul + add_assign", move |b| {
        b.iter(|| mul_then_add(black_box(&bases_clone), black_box(&scalars_clone)))
    });
}

fn criterion_benchmark(c: &mut Criterion) {
    msm::<ark_bw6_761::G1Affine>(c, 6);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);