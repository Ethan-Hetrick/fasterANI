//! Runtime options plus memory/RSS accounting and progress-reporting helpers.

#[cfg(debug_assertions)]
use std::{env, mem::size_of};
use std::{fs, time::Instant};

use crate::ani::CliArgs;
#[cfg(debug_assertions)]
use crate::ani::{MinimizerKey, ReferenceMinimizer, SeedHit};

/// Runtime controls shared by long-running reference build and sketch operations.
#[derive(Clone, Copy, Default)]
pub(crate) struct RuntimeOptions {
    pub(crate) progress_enabled: bool,
    pub(crate) worker_threads: usize,
}

impl RuntimeOptions {
    pub(crate) fn effective_worker_threads(self) -> usize {
        self.worker_threads.max(1)
    }

    pub(crate) fn with_worker_threads(self, worker_threads: usize) -> Self {
        Self {
            worker_threads: worker_threads.max(1),
            ..self
        }
    }
}

#[cfg(debug_assertions)]
#[derive(Default)]
pub(crate) struct QueryMemoryEstimate {
    pub(crate) fragment_struct_bytes: usize,
    pub(crate) query_minimizer_vec_bytes: usize,
    pub(crate) seed_minimizer_vec_bytes: usize,
}

/// Return the current process peak RSS in kilobytes, or -1 when unavailable.
pub(crate) fn peak_rss_kb() -> i64 {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    let result = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
    if result == 0 {
        unsafe { usage.assume_init().ru_maxrss }
    } else {
        -1
    }
}

pub(crate) fn memory_mib(bytes: usize) -> f64 {
    bytes as f64 / 1024.0 / 1024.0
}

fn current_rss_kb() -> i64 {
    let Ok(status) = fs::read_to_string("/proc/self/status") else {
        return peak_rss_kb();
    };

    for line in status.lines() {
        let Some(rest) = line.strip_prefix("VmRSS:") else {
            continue;
        };
        let Some(value) = rest.split_whitespace().next() else {
            continue;
        };
        if let Ok(kb) = value.parse::<i64>() {
            return kb;
        }
    }

    peak_rss_kb()
}

pub(crate) fn emit_progress(stage: &str, message: &str, start: Instant) {
    let rss_kb: i64 = current_rss_kb();
    let rss_gib: f64 = if rss_kb > 0 {
        rss_kb as f64 / 1024.0 / 1024.0
    } else {
        f64::NAN
    };
    eprintln!(
        "PROGRESS\tstage={stage}\t{message}\trss_gib={rss_gib:.3}\telapsed_s={:.3}",
        start.elapsed().as_secs_f64()
    );
}

#[cfg(debug_assertions)]
pub(crate) fn reference_build_struct_bytes(
    reference_minimizers: usize,
    seed_hits: usize,
    unique_index_keys: usize,
) -> usize {
    reference_minimizers
        .saturating_mul(size_of::<ReferenceMinimizer>())
        .saturating_add(seed_hits.saturating_mul(size_of::<SeedHit>()))
        .saturating_add(unique_index_keys.saturating_mul(size_of::<MinimizerKey>()))
}

#[cfg(debug_assertions)]
pub(crate) fn performance_metrics_enabled(args: &CliArgs) -> bool {
    args.verbose || env::var_os("FASTERANI_METRICS").is_some()
}

#[cfg(not(debug_assertions))]
pub(crate) fn performance_metrics_enabled(_args: &CliArgs) -> bool {
    false
}
