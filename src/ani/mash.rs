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

pub(crate) fn binomial_survival(x: usize, probability: f64, trials: usize) -> f64 {
    binomial_survival_at_least(x, probability, trials)
}

/// Approximate t-distribution CDF. Returns P(X <= t) for X ~ t(df).
pub(crate) fn t_cdf_approx(t: f64, df: f64) -> f64 {
    if !t.is_finite() || !df.is_finite() || df <= 0.0 {
        return f64::NAN;
    }

    if df > 30.0 {
        return standard_normal_cdf(t);
    }

    let x: f64 = df / (df + t * t);
    let beta: f64 = regularized_incomplete_beta(x, df / 2.0, 0.5);

    if t >= 0.0 {
        1.0 - 0.5 * beta
    } else {
        0.5 * beta
    }
    .clamp(0.0, 1.0)
}

fn standard_normal_cdf(value: f64) -> f64 {
    (0.5 * (1.0 + erf(value / 2.0_f64.sqrt()))).clamp(0.0, 1.0)
}

fn erf(value: f64) -> f64 {
    let sign: f64 = if value < 0.0 { -1.0 } else { 1.0 };
    let x: f64 = value.abs();
    let t: f64 = 1.0 / (1.0 + 0.327_591_1 * x);
    let y: f64 = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t
            + 0.254_829_592)
            * t
            * (-x * x).exp();

    sign * y
}

fn regularized_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }

    let bt: f64 =
        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();

    if x < (a + 1.0) / (a + b + 2.0) {
        bt * beta_continued_fraction(x, a, b) / a
    } else {
        1.0 - bt * beta_continued_fraction(1.0 - x, b, a) / b
    }
}

fn beta_continued_fraction(x: f64, a: f64, b: f64) -> f64 {
    const MAX_ITERATIONS: usize = 200;
    const EPSILON: f64 = 3.0e-14;
    const MIN_FLOAT: f64 = 1.0e-300;

    let qab: f64 = a + b;
    let qap: f64 = a + 1.0;
    let qam: f64 = a - 1.0;
    let mut c: f64 = 1.0;
    let mut d: f64 = 1.0 - qab * x / qap;
    if d.abs() < MIN_FLOAT {
        d = MIN_FLOAT;
    }
    d = 1.0 / d;
    let mut h: f64 = d;

    for m in 1..=MAX_ITERATIONS {
        let m_f: f64 = m as f64;
        let m2: f64 = 2.0 * m_f;

        let mut aa: f64 = m_f * (b - m_f) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < MIN_FLOAT {
            d = MIN_FLOAT;
        }
        c = 1.0 + aa / c;
        if c.abs() < MIN_FLOAT {
            c = MIN_FLOAT;
        }
        d = 1.0 / d;
        h *= d * c;

        aa = -(a + m_f) * (qab + m_f) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < MIN_FLOAT {
            d = MIN_FLOAT;
        }
        c = 1.0 + aa / c;
        if c.abs() < MIN_FLOAT {
            c = MIN_FLOAT;
        }
        d = 1.0 / d;
        let delta: f64 = d * c;
        h *= delta;

        if (delta - 1.0).abs() < EPSILON {
            break;
        }
    }

    h
}

fn ln_gamma(value: f64) -> f64 {
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];

    if value < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * value).sin().ln()
            - ln_gamma(1.0 - value);
    }

    let z: f64 = value - 1.0;
    let mut x: f64 = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        x += coefficient / (z + index as f64);
    }
    let t: f64 = z + 7.5;

    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
}

pub(crate) fn mash_distance_lower_bound(
    distance: f64,
    sketch_size: usize,
    kmer_size: usize,
    confidence_interval: f64,
) -> f64 {
    if confidence_interval <= 0.0 {
        return distance;
    }

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
    x = x.min(sketch_size);

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
    if percent_identity <= 0.0 {
        return 1;
    }

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

#[cfg(test)]
mod tests {
    use super::{mash_distance_lower_bound, t_cdf_approx};

    #[test]
    fn t_cdf_approx_matches_table_value() {
        assert!((t_cdf_approx(2.262, 9.0) - 0.975).abs() < 0.01);
    }

    #[test]
    fn mash_distance_lower_bound_handles_confidence_endpoints() {
        let distance = 0.1;
        assert_eq!(mash_distance_lower_bound(distance, 1_000, 16, 0.0), distance);

        let full_confidence_bound = mash_distance_lower_bound(distance, 1_000, 16, 1.0);
        assert!(full_confidence_bound.is_finite());
        assert!(full_confidence_bound <= distance);
    }
}
