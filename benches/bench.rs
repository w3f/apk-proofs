use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ark_ff::{Field, PrimeField, Zero, test_rng, UniformRand};
use ark_ec::{AffineCurve, ProjectiveCurve};

extern crate apk_proofs;

fn msm<G: AffineCurve>(c: &mut Criterion, n: usize) {

    let rng = &mut test_rng();

    let nu = G::ScalarField::rand(rng);
    let scalars = (0..n).map(|i| nu.pow([i as u64]).into_repr()).collect::<Vec<_>>();
    let bases = (0..n).map(|_| G::Projective::rand(rng).into_affine()).collect::<Vec<_>>();

    {
        let (scalars, bases) = (scalars.clone(), bases.clone());
        c.bench_function("ark_ec::msm::VariableBaseMSM", move |b| {
            b.iter(|| ark_ec::msm::VariableBaseMSM::multi_scalar_mul(black_box(&bases), black_box(&scalars)))
        });
    }

    {
        let (scalars, bases) = (scalars.clone(), bases.clone());
        c.bench_function("naive mul then add", move |b| {
            b.iter(|| apk_proofs::utils::mul_then_add(black_box(&bases), black_box(&scalars)))
        });
    }

    let nu: u128 = rand::random();
    let nu = G::ScalarField::from(nu);

    {
        let bases = bases.clone();
        let nu = nu.into_repr();

        c.bench_function("128-bit Horner", move |b| {
            b.iter(|| apk_proofs::utils::horner(black_box(&bases), black_box(nu)))
        });
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    msm::<ark_bw6_761::G1Affine>(c, 6);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);