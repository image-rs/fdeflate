#![feature(test)]

extern crate test;

use fdeflate::compress_to_vec;
use rand::RngExt;

#[bench]
fn bench_compute_code_lengths(b: &mut test::Bencher) {
    const N: usize = 48;
    let mut rng = rand::rng();
    let mut freqs = vec![0; N];
    for f in freqs.iter_mut() {
        *f = rng.random_range::<u64, _>(1..1000);
    }

    b.iter(|| {
        let mut lengths = vec![0; N];
        fdeflate::compute_code_lengths(&freqs, &[1; N], &[8; N], &mut lengths);
    });
}

#[bench]
fn bench_uniform_random(b: &mut test::Bencher) {
    let mut rng = rand::rng();
    let mut data = vec![0; 1024 * 1024];
    for byte in &mut data {
        *byte = rng.random();
    }
    b.bytes = data.len() as u64;
    b.iter(|| compress_to_vec(&data));
}

#[bench]
fn bench_low(b: &mut test::Bencher) {
    let mut rng = rand::rng();
    let mut data = vec![0; 1024 * 1024];
    for byte in &mut data {
        *byte = (rng.random_range::<u8, _>(0..16) * 2).wrapping_sub(16);
    }
    b.bytes = data.len() as u64;
    b.iter(|| compress_to_vec(&data));
}

#[bench]
fn bench_mixture(b: &mut test::Bencher) {
    let mut rng = rand::rng();
    let mut data = vec![0; 1024 * 1024];
    for byte in &mut data {
        if rng.random_range(0..200) == 1 {
            *byte = rng.random();
        } else {
            *byte = rng.random_range::<u8, _>(0..32).wrapping_sub(16);
        }
    }
    b.bytes = data.len() as u64;
    b.iter(|| compress_to_vec(&data));
}

#[bench]
fn bench_distribution(b: &mut test::Bencher) {
    let mut rng = rand::rng();
    let mut data = vec![0; 1024 * 1024];
    for byte in &mut data {
        *byte = match rng.random_range(0..100) {
            0 => rng.random(),
            1..=2 => rng.random_range::<u8, _>(0..32).wrapping_sub(16),
            11..=50 => rng.random_range::<u8, _>(0..16).wrapping_sub(8),
            51..=80 => rng.random_range::<u8, _>(0..8).wrapping_sub(4),
            _ => 0,
        }
    }
    b.bytes = data.len() as u64;
    b.iter(|| compress_to_vec(&data));
}
