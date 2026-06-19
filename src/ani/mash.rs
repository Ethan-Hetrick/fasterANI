//! Mash-distance math and shared-minimizer lower-bound estimates.

pub(crate) fn fastani_mash_distance(jaccard: f64, kmer_size: usize) -> f64 {
    if jaccard == 0.0 {
        return 1.0;
    }

    if jaccard == 1.0 {
        return 0.0;
    }

    (-1.0 / kmer_size as f64) * ((2.0 * jaccard) / (1.0 + jaccard)).ln()
}

fn mash_distance_to_jaccard(distance: f64, kmer_size: usize) -> f64 {
    1.0 / (2.0 * (kmer_size as f64 * distance).exp() - 1.0)
}

fn binomial_survival_at_least(x: usize, probability: f64, trials: usize) -> f64 {
    if x == 0 {
        return 1.0;
    }

    if x > trials || probability <= 0.0 {
        return 0.0;
    }

    if probability >= 1.0 {
        return 1.0;
    }

    let q = 1.0 - probability;
    let mut pmf = q.powi(trials as i32);
    let mut cdf_below = pmf;

    for i in 0..(x - 1) {
        pmf *= (trials - i) as f64 / (i + 1) as f64 * probability / q;
        cdf_below += pmf;
    }

    (1.0 - cdf_below).clamp(0.0, 1.0)
}

pub(crate) fn mash_distance_lower_bound(
    distance: f64,
    sketch_size: usize,
    kmer_size: usize,
    confidence_interval: f64,
) -> f64 {
    let q2: f64 = (1.0 - confidence_interval) / 2.0;
    let jaccard: f64 = mash_distance_to_jaccard(distance, kmer_size);
    let mut x: usize = ((sketch_size as f64 * jaccard).ceil() as usize).max(1);

    while x <= sketch_size {
        let cdf_complement = binomial_survival_at_least(x, jaccard, sketch_size);

        if cdf_complement < q2 {
            x = x.saturating_sub(1);
            break;
        }

        x += 1;
    }

    fastani_mash_distance(x as f64 / sketch_size as f64, kmer_size)
}

fn estimate_minimum_shared_minimizers(
    sketch_size: usize,
    kmer_size: usize,
    percent_identity: f64,
) -> usize {
    let mash_distance: f64 = 1.0 - percent_identity / 100.0;
    let jaccard: f64 = mash_distance_to_jaccard(mash_distance, kmer_size);

    (sketch_size as f64 * jaccard).ceil() as usize
}

/// Estimate the relaxed minimum shared minimizers needed to keep a candidate alive.
pub(crate) fn estimate_relaxed_minimum_shared_minimizers(
    sketch_size: usize,
    kmer_size: usize,
    percent_identity: f64,
    mash_confidence: f64,
) -> usize {
    let strict_minimum: usize =
        estimate_minimum_shared_minimizers(sketch_size, kmer_size, percent_identity);
    let mut relaxed_minimum: usize = strict_minimum;

    for i in (0..=strict_minimum).rev() {
        let jaccard = i as f64 / sketch_size as f64;
        let distance = fastani_mash_distance(jaccard, kmer_size);
        let lower_distance =
            mash_distance_lower_bound(distance, sketch_size, kmer_size, mash_confidence);
        let upper_identity = 100.0 * (1.0 - lower_distance);

        if upper_identity >= percent_identity {
            relaxed_minimum = i;
        } else {
            break;
        }
    }

    relaxed_minimum
}
