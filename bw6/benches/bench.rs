use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use ark_ff::{Field, PrimeField};
use ark_std::{UniformRand, test_rng};
use ark_ec::{AffineCurve, ProjectiveCurve};

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

fn verification(c: &mut Criterion) {
    use apk_proofs::{Prover, Verifier, Setup, SignerSet, Bitmask, bls};
    use merlin::Transcript;
    use rand::{Rng, seq::SliceRandom};
    use std::convert::TryInto;

    let mut group = c.benchmark_group("verification");

    let rng = &mut test_rng();
    let log_domain_size_range = 8..=10;

    for log_domain_size in log_domain_size_range {
        let setup = Setup::generate(log_domain_size, rng);
        let keyset_size = rng.gen_range(1..=setup.max_keyset_size());
        let keyset_size = keyset_size.try_into().unwrap();
        let signer_set = SignerSet::random(keyset_size, rng);
        let pks_comm = signer_set.commit(setup.domain_size, &setup.kzg_params.get_pk());
        let bits = (0..keyset_size).map(|_| rng.gen_bool(2.0 / 3.0)).collect::<Vec<_>>();
        let bitmask = Bitmask::from_bits(&bits);
        let apk = bls::PublicKey::aggregate(signer_set.get_by_mask(&bitmask));

        let prover = Prover::new(
            setup.domain_size,
            setup.kzg_params.get_pk(),
            &pks_comm,
            signer_set.get_all(),
            Transcript::new(b"apk_proof"),
        );
        let proof_basic = prover.prove_simple(bitmask.clone());
        let proof_packed = prover.prove_packed(bitmask.clone());
        let proof_counting = prover.prove_counting(bitmask.clone());

        let create_verifier = || {
            Verifier::new(
                setup.domain_size,
                setup.kzg_params.get_vk(),
                pks_comm.clone(),
                Transcript::new(b"apk_proof"),
            )
        };

        group.bench_with_input(
            BenchmarkId::new("basic", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                let verifier = create_verifier();
                verifier.verify_simple(&apk, &bitmask, &proof_basic);
            }),
        );

        group.bench_with_input(
            BenchmarkId::new("packed", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                let verifier = create_verifier();
                verifier.verify_packed(&apk, &bitmask, &proof_packed);
            }),
        );

        group.bench_with_input(
            BenchmarkId::new("counting", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| {
                let verifier = create_verifier();
                verifier.verify_counting(&apk, &bitmask, &proof_counting);
            }),
        );
    }
    group.finish();
}

fn components(c: &mut Criterion) {
    msm::<ark_bw6_761::G1Affine>(c, 6);
    barycentric_evaluation::<ark_bw6_761::Fr>(c, 2u32.pow(10));
    bw6_subgroup_check(c);
}

criterion_group!(benches, components, verification);
criterion_main!(benches);