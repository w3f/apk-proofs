use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ark_ff::{Field, PrimeField, test_rng, UniformRand};
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

fn apk_verification(c: &mut Criterion) {
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain};
    use apk_proofs::Params;
    use bitvec::vec::BitVec;
    use rand::Rng;

    let num_pks = 1000;

    let rng = &mut test_rng();

    let signer_set = apk_proofs::SignerSet::random(num_pks, rng);
    let params = Params::new(signer_set.size(), rng);
    let pks_domain_size = GeneralEvaluationDomain::<ark_bw6_761::Fr>::compute_size_of_domain(num_pks).unwrap();
    let (pks_x_comm, pks_y_comm) = signer_set.commit(&params.get_ck(pks_domain_size));
    let bitmask: BitVec = (0..num_pks).map(|_| rng.gen_bool(2.0 / 3.0)).collect();
    let apk = apk_proofs::bls::PublicKey::aggregate(signer_set.get_by_mask(&bitmask));
    let proof = apk_proofs::prove(&bitmask, signer_set.get_all(), &params.to_pk(), ProofScheme::Accountable);
    let vk = params.to_vk();
    c.bench_function("apk verification", move |b| {
        b.iter(|| apk_proofs::verify(
            black_box(&pks_x_comm),
            black_box(&pks_y_comm),
            black_box(&apk),
            black_box(&bitmask),
            black_box(&proof),
            black_box(&vk)
        ))
    });
}

fn criterion_benchmark(c: &mut Criterion) {
    msm::<ark_bw6_761::G1Affine>(c, 6);
    barycentric_evaluation::<ark_bw6_761::Fr>(c, 2u32.pow(10));
    bw6_subgroup_check(c);
    apk_verification(c);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
