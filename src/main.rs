use std::fs;
use std::io::Write;

fn main() {
    let mut total_bytes = 0;
    let mut compressed_bytes = 0;
    let mut ratios = Vec::new();

    let mut zero_run = 0u64;
    let mut run = 0u64;
    let mut pixel_run = 0u64;

    let mut counts = [0u64; 286];
    let mut total_counts = [1u64; 286];

    for (i, entry) in fs::read_dir("/home/jonathan/git/image-png/rawdata2")
        .unwrap()
        .flatten()
        .enumerate()
    {
        let raw = fs::read(entry.path()).unwrap();
        let data = fdeflate::compress_to_vec(&raw);
        let ratio = data.len() as f64 / raw.len() as f64;
        ratios.push(ratio * 100.0);

        total_bytes += raw.len();
        compressed_bytes += data.len();

        if i % 10 == 9 {
            print!(".");
            std::io::stdout().flush().unwrap();
            if i % 500 == 499 {
                println!();
            }
        }

        let mut current_run = 0;

        let mut i = 0;
        while i < raw.len() || current_run > 0 {
            if raw.get(i) == Some(&0) {
                current_run += 1;
            } else {
                // Runs start with a literal
                if current_run > 0 {
                    counts[0] += 1;
                    current_run -= 1;
                }

                // Run component
                while current_run >= 258 {
                    counts[285] += 1;
                    current_run -= 258;
                }
                if current_run >= 5 {
                    counts[fdeflate::LEN_SYM[current_run as usize - 3] as usize] += 1;
                } else {
                    counts[0] += current_run as u64;
                }
                current_run = 0;

                // Literal of the symbol that ended the run
                if let Some(&b) = raw.get(i) {
                    counts[b as usize] += 1;
                }
            }
            i += 1;
        }

        for i in 3..raw.len() {
            if raw[i - 1] == raw[i] && raw[i - 2] == raw[i] && raw[i - 3] == raw[i] {
                if raw[i] == 0 {
                    zero_run += 1;
                } else {
                    run += 1;
                }
            }
        }

        for a in raw.chunks_exact(9) {
            if a != [0; 9] && a[0..6] == a[3..9] {
                pixel_run += 9;
            }
        }

        for i in 0..286 {
            total_counts[i] += counts[i] /* * (1<<32) / raw.len() as u64*/;
        }
        counts = [0; 286];
    }
    println!();
    println!();

    counts = total_counts;

    let mut lengths = vec![0; 286];
    fdeflate::compute_code_lengths(&counts, &[1; 286], &[12; 286], &mut lengths);
    let compression12 = lengths
        .iter()
        .zip(&counts)
        .map(|(&l, &c)| u64::from(l) * c)
        .sum::<u64>() as f64
        / 8.0
        / total_bytes as f64;

    println!();
    println!("lengths = {:?}", lengths);
    println!();

    fdeflate::compute_code_lengths(&counts, &[1; 286], &[15; 286], &mut lengths);
    let compression15 = lengths
        .iter()
        .zip(&counts)
        .map(|(&l, &c)| u64::from(l) * c)
        .sum::<u64>() as f64
        / 8.0
        / total_bytes as f64;

    println!(
        "compression12 = {}%, compression15 = {}%",
        compression12 * 100.0,
        compression15 * 100.0
    );
    println!();

    println!(
        "Zero run: {:.2}%, Nonzero: {:.2}%",
        100.0 * zero_run as f64 / total_bytes as f64,
        100.0 * run as f64 / total_bytes as f64
    );
    println!(
        "Pixel run: {:.2}%",
        100.0 * pixel_run as f64 / total_bytes as f64
    );

    let mean = ratios.iter().sum::<f64>() / ratios.len() as f64;
    let geomean = (ratios.iter().map(|r| r.ln()).sum::<f64>() / ratios.len() as f64).exp();
    println!(
        "total: {:.2}%, mean: {:.2}%, geomean: {:.2}%",
        100.0 * compressed_bytes as f64 / total_bytes as f64,
        mean,
        geomean
    );
}
