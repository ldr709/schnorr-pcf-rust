use core::time::Duration;
use criterion::*;
use pcf_threshold_schnorr::{dl_pcf::*, schnorr_2p::*};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[allow(non_snake_case)]
pub fn benchmark_dl_pcf(c: &mut Criterion) {
    let mut rng = ChaCha20Rng::from_entropy();

    // Disable the very slow benchmarks for now.
    if false {
        let mut group = c.benchmark_group("Setup functions");
        group.sampling_mode(SamplingMode::Flat);

        // Minimum settings for iteration count in Criterion.
        group.sample_size(10);
        group.nresamples(10);
        group.measurement_time(Duration::from_secs(1));

        group.bench_function("DL PCF Setup", |b| b.iter_with_large_drop(|| setup(&mut rng, 3072)));
        group.bench_function("2-party Schorr keygen", |b| b.iter_with_large_drop(|| keygen(&mut rng)));
        group.finish();
    }

    let (prover, verifier) = setup(&mut rng, 3072);
    let msg = b"Benchmarking...";
    c.bench_function("DL PCF Prove", |b| b.iter_with_large_drop(||
        prover.prove(msg)));

    let (_r, _R, proof) = prover.prove(msg);
    c.bench_function("DL PCF Verify", |b| b.iter_with_large_drop(||
        verifier.verify(msg, proof.clone())));

    let (_pk, signers) = keygen(&mut rng);
    c.bench_function("2-party Schnorr", |b| b.iter_with_large_drop(||
        sign_2party(&signers, msg)));
}

criterion_group!(benches, benchmark_dl_pcf);
criterion_main!(benches);
