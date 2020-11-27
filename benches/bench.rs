use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ark_ff::{Field, PrimeField, test_rng, UniformRand};
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

    let nu: G::ScalarField = u128::rand(rng).into();

    {
        let bases = bases.clone();
        let nu = nu.into_repr();

        c.bench_function("128-bit Horner", move |b| {
            b.iter(|| apk_proofs::utils::horner(black_box(&bases), black_box(nu)))
        });
    }
}

fn bw6_subgroup_check(c: &mut Criterion) {
    let rng = &mut test_rng();

    let p = ark_bw6_761::G1Projective::rand(rng);
    let p_affine = p.into_affine();

    c.bench_function("subgroup check: mul by group order", move |b| {
        b.iter(|| (black_box(&p_affine)).is_in_correct_subgroup_assuming_on_curve())
    });

    c.bench_function("subgroup check: GLV", move |b| {
        b.iter(|| apk_proofs::endo::subgroup_check(black_box(&p)))
    });
}

fn criterion_benchmark(c: &mut Criterion) {
    msm::<ark_bw6_761::G1Affine>(c, 6);
    bw6_subgroup_check(c);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);