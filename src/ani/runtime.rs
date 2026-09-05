//! Runtime options plus memory/RSS accounting and progress-reporting helpers.

#[cfg(debug_assertions)]
use std::{env, mem::size_of};
use std::{fs, io, time::Instant};

use crate::ani::{cli::CliArgs, constants::DEFAULT_MPHF_GAMMA};
#[cfg(debug_assertions)]
use crate::ani::{
    constants::MinimizerKey,
    model::reference::{ReferenceMinimizer, SeedHit},
};

/// Runtime controls shared by long-running reference build and sketch operations.
#[derive(Clone, Copy)]
pub(crate) struct RuntimeOptions {
    pub(crate) progress_enabled: bool,
    pub(crate) worker_threads: usize,
    pub(crate) mphf_gamma: f64,
    build_progress: Option<BuildProgressContext>,
}

#[derive(Clone, Copy)]
struct BuildProgressContext {
    generation_id: [u8; 64],
    generation_id_len: u8,
    shard_index: usize,
}

impl Default for RuntimeOptions {
    fn default() -> Self {
        Self {
            progress_enabled: false,
            worker_threads: 0,
            mphf_gamma: DEFAULT_MPHF_GAMMA,
            build_progress: None,
        }
    }
}

impl RuntimeOptions {
    pub(crate) fn with_progress_enabled(self, progress_enabled: bool) -> Self {
        Self {
            progress_enabled,
            ..self
        }
    }

    pub(crate) fn effective_worker_threads(self) -> usize {
        self.worker_threads.max(1)
    }

    pub(crate) fn with_worker_threads(self, worker_threads: usize) -> Self {
        Self {
            worker_threads: worker_threads.max(1),
            ..self
        }
    }

    pub(crate) fn with_mphf_gamma(self, mphf_gamma: f64) -> Self {
        Self { mphf_gamma, ..self }
    }

    pub(crate) fn with_build_progress(
        self,
        generation_id: &str,
        shard_index: usize,
    ) -> io::Result<Self> {
        let generation_id_len: u8 = u8::try_from(generation_id.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "database build generation identifier is too long for progress attribution",
            )
        })?;
        let mut generation_id_bytes: [u8; 64] = [0; 64];
        let destination = generation_id_bytes
            .get_mut(..generation_id.len())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "database build generation identifier exceeds 64 bytes",
                )
            })?;
        destination.copy_from_slice(generation_id.as_bytes());
        Ok(Self {
            build_progress: Some(BuildProgressContext {
                generation_id: generation_id_bytes,
                generation_id_len,
                shard_index,
            }),
            ..self
        })
    }
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

/// Emit a nested build event with explicit generation, shard, and assigned-worker attribution.
pub(crate) fn emit_runtime_progress(
    runtime_options: RuntimeOptions,
    stage: &str,
    message: &str,
    start: Instant,
) {
    if let Some(context) = runtime_options.build_progress {
        let generation_id =
            std::str::from_utf8(&context.generation_id[..usize::from(context.generation_id_len)])
                .expect("build progress generation id came from UTF-8");
        emit_progress(
            stage,
            &format!(
                "generation_id={generation_id}\tshard={}\tassigned_threads={}\t{message}",
                context.shard_index,
                runtime_options.effective_worker_threads()
            ),
            start,
        );
    } else {
        emit_progress(stage, message, start);
    }
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
