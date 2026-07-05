use std::{process::ExitCode, time::Instant};

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

fn main() -> ExitCode {
    let total_start: Instant = Instant::now();

    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    // Print errors with Display (real newlines, no `Custom { .. }` debug noise) and
    // exit non-zero, rather than letting the runtime Debug-format the io::Error.
    match faster_ani::run_started_at(total_start) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
