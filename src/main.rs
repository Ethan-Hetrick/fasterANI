use std::process::ExitCode;

fn main() -> ExitCode {
    // Print errors with Display (real newlines, no `Custom { .. }` debug noise) and
    // exit non-zero, rather than letting the runtime Debug-format the io::Error.
    match faster_ani::run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
