use criterion::{black_box, criterion_group, criterion_main, Criterion};

extern crate apk_proofs;

fn criterion_benchmark(c: &mut Criterion) {

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
criterion_main!(benches);