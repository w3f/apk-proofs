use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use ark_ff::{Field, PrimeField, One, FftField};
use ark_std::{UniformRand, test_rng};
use ark_ec::{AffineCurve, ProjectiveCurve};
use apk_proofs::{Keyset, setup};
use ark_bw6_761::Fr;
use ark_poly::{Evaluations, EvaluationDomain, UVPolynomial, Radix2EvaluationDomain};
use ark_poly::univariate::DensePolynomial;

extern crate apk_proofs;

fn barycentric_evaluation<F: Field>(c: &mut Criterion, n: u32) {
    use ark_poly::{Evaluations, EvaluationDomain, Radix2EvaluationDomain, Polynomial};

    let rng = &mut test_rng();
    let n = std::convert::TryInto::try_into(n).unwrap();
    let domain = Radix2EvaluationDomain::new(n).unwrap();
    let evals = (0..n).map(|_| ark_bw6_761::Fr::rand(rng)).collect::<Vec<_>>();
    let evals2 = evals.clone();
    let z = ark_bw6_761::Fr::rand(rng);

    c.bench_function("barycentric_evaluation", move |b| {
        b.iter(|| {
            apk_proofs::utils::barycentric_eval_at(black_box(z), black_box(&evals), black_box(domain))
        })
    });

    let evals = Evaluations::from_vec_and_domain(evals2, domain);
    c.bench_function("interpolate + evaluate", move |b| {
        b.iter(|| {
            black_box(&evals).interpolate_by_ref().evaluate(black_box(&z));
        })
    });
}


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

fn amplification(c: &mut Criterion) {
    use ark_poly::{EvaluationDomain, Radix2EvaluationDomain, UVPolynomial};
    use apk_proofs::domains::Domains;

    let mut group = c.benchmark_group("amplification");

    let rng = &mut test_rng();
    let log_domain_size_range = vec![10, 16];

    for log_domain_size in log_domain_size_range {
        let n = 2u32.pow(log_domain_size) as usize;
        let domains = Domains::new(n);

        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();

        group.bench_with_input(
            BenchmarkId::new("2x", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                let poly = Evaluations::from_vec_and_domain(evals.clone(), domains.domain).interpolate();
                poly.evaluate_over_domain_by_ref(domains.domain2x)
            }),
        );

        group.bench_with_input(
            BenchmarkId::new("2x-coset", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                domains.amplify_x2(evals.clone())
            }),
        );

        group.bench_with_input(
            BenchmarkId::new("4x", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                let poly = Evaluations::from_vec_and_domain(evals.clone(), domains.domain).interpolate();
                poly.evaluate_over_domain_by_ref(domains.domain4x)
            }),
        );

        group.bench_with_input(
            BenchmarkId::new("4x-coset", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                domains.amplify_x4(evals.clone())
            }),
        );
    }
}

fn verification(c: &mut Criterion) {
    use apk_proofs::{Prover, Verifier, Bitmask, bls};
    use merlin::Transcript;
    use rand::{Rng, seq::SliceRandom};
    use std::convert::TryInto;

    let mut group = c.benchmark_group("verification");

    let rng = &mut test_rng();
    let log_domain_size_range = 8..=10;

    for log_domain_size in log_domain_size_range {
        let keyset_size = (2u32.pow(log_domain_size) - 1) as usize;
        let keyset = (0..keyset_size).map(|_| ark_bls12_377::G1Projective::rand(rng)).collect();
        let keyset = Keyset::new(keyset);
        let kzg_params = setup::generate_for_keyset(keyset_size, rng);
        let pks_comm = keyset.commit(&kzg_params.get_pk());

        let bitmask = Bitmask::from_bits(&vec![true; keyset_size]);

        let prover = Prover::new(
            keyset,
            &pks_comm,
            kzg_params.clone(),
            Transcript::new(b"apk_proof"),
        );
        let proof_basic = prover.prove_simple(bitmask.clone());
        let proof_packed = prover.prove_packed(bitmask.clone());
        let proof_counting = prover.prove_counting(bitmask.clone());

        let create_verifier = || {
            Verifier::new(
                kzg_params.get_vk(),
                pks_comm.clone(),
                Transcript::new(b"apk_proof"),
            )
        };

        group.bench_with_input(
            BenchmarkId::new("basic", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                let verifier = create_verifier();
                verifier.verify_simple(&proof_basic.1, &proof_basic.0);
            }),
        );

        group.bench_with_input(
            BenchmarkId::new("packed", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                let verifier = create_verifier();
                verifier.verify_packed(&proof_packed.1, &proof_packed.0);
            }),
        );

        let count = bitmask.count_ones();
        group.bench_with_input(
            BenchmarkId::new("counting", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                let verifier = create_verifier();
                verifier.verify_counting(&proof_counting.1, &proof_counting.0);
            }),
        );
    }

    group.finish();
}

fn fft<F: FftField, D: EvaluationDomain<F>>(c: &mut Criterion) {
    let mut group = c.benchmark_group("FFT");

    let rng = &mut test_rng();

    for logn in 10..=16 {
        let n = 2usize.pow(logn);
        let domain = D::new(n).unwrap();
        let poly = DensePolynomial::<F>::rand(n-1, rng);
        let coeffs = poly.coeffs;

        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, n| {
            b.iter(|| domain.fft(&coeffs))
        });
    }
    group.finish();
}


fn components(c: &mut Criterion) {
    msm::<ark_bw6_761::G1Affine>(c, 6);
    barycentric_evaluation::<ark_bw6_761::Fr>(c, 2u32.pow(10));
    bw6_subgroup_check(c);
    amplification(c);
}

fn primitives(c: &mut Criterion) {
    fft::<ark_bw6_761::Fr, Radix2EvaluationDomain<ark_bw6_761::Fr>>(c);
}

criterion_group!(benches, components, verification, primitives);
criterion_main!(benches);