//! Loading and validating TOML runtime parameter files.

use std::{fs, io};

use serde::Deserialize;

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ParamsFileConfig {
    pub(crate) threads: Option<usize>,
    pub(crate) freq_threshold_percent: Option<f64>,
    pub(crate) minmer_count: Option<usize>,
    pub(crate) kmer_size: Option<usize>,
    pub(crate) window_size: Option<usize>,
    pub(crate) fragment_length: Option<u32>,
    pub(crate) fragment_stride: Option<u32>,
    pub(crate) min_fragment_length: Option<u32>,
    pub(crate) mash_threshold: Option<f64>,
    pub(crate) mash_confidence: Option<f64>,
    pub(crate) mphf_gamma: Option<f64>,
    pub(crate) split_n_run: Option<usize>,
    pub(crate) max_shard_minimizers: Option<usize>,
    pub(crate) shards: Option<String>,
    pub(crate) index_build_mode: Option<String>,
    pub(crate) bgzip: Option<bool>,
    pub(crate) header: Option<bool>,
    pub(crate) verbose: Option<bool>,
    pub(crate) quiet: Option<bool>,
    pub(crate) force: Option<bool>,
    pub(crate) reference_files: Option<Vec<String>>,
    pub(crate) reference_lists: Option<Vec<String>>,
    pub(crate) query_files: Option<Vec<String>>,
    pub(crate) query_lists: Option<Vec<String>>,
    pub(crate) query_name: Option<String>,
    pub(crate) reference_sketch: Option<String>,
    pub(crate) tmp: Option<String>,
    pub(crate) out: Option<String>,
    pub(crate) mapping_stats: Option<String>,
}

pub(crate) fn load_params_file(path: &str) -> io::Result<ParamsFileConfig> {
    validate_params_file_path(path)?;
    let contents = fs::read_to_string(path).map_err(|err| {
        io::Error::new(
            err.kind(),
            format!("Failed to read params file at {path}: {err}"),
        )
    })?;

    toml::from_str::<ParamsFileConfig>(&contents).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Failed to parse params file at {path}:\n{err}"),
        )
    })
}

pub(crate) fn validate_params_file_path(path: &str) -> io::Result<()> {
    let metadata = fs::metadata(path).map_err(|err| {
        io::Error::new(
            err.kind(),
            format!("Params file not found or unreadable: {path}\n{err}"),
        )
    })?;
    if !metadata.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Params file path is not a file: {path}"),
        ));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        time::{SystemTime, UNIX_EPOCH},
    };

    use super::{load_params_file, ParamsFileConfig};

    fn temp_params_path(name: &str) -> String {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before unix epoch")
            .as_nanos();
        std::env::temp_dir()
            .join(format!("fasterani-{name}-{nanos}.toml"))
            .to_string_lossy()
            .into_owned()
    }

    #[test]
    fn valid_toml_file_loads_successfully() {
        let path = temp_params_path("valid");
        fs::write(
            &path,
            r#"
threads = 4
mash_threshold = 82.5
mash_confidence = 0.79
mphf_gamma = 2.0
reference_files = ["ref.fa"]
"#,
        )
        .expect("write temp params file");

        let config = load_params_file(&path).expect("load params file");
        assert_eq!(config.threads, Some(4));
        assert_eq!(config.mash_threshold, Some(82.5));
        assert_eq!(config.mash_confidence, Some(0.79));
        assert_eq!(config.mphf_gamma, Some(2.0));
        assert_eq!(config.reference_files, Some(vec!["ref.fa".to_owned()]));

        let _ = fs::remove_file(path);
    }

    #[test]
    fn malformed_toml_returns_clear_error() {
        let path = temp_params_path("malformed");
        fs::write(&path, "threads = true42").expect("write temp params file");

        let err = load_params_file(&path).expect_err("malformed file should fail");
        assert!(err.to_string().contains("Failed to parse params file"));

        let _ = fs::remove_file(path);
    }

    #[test]
    fn missing_file_returns_file_not_found() {
        let path = temp_params_path("missing");
        let err = load_params_file(&path).expect_err("missing file should fail");
        assert!(err.to_string().contains("Params file not found"));
    }

    #[test]
    fn optional_fields_default_to_none() {
        let config: ParamsFileConfig = toml::from_str("").expect("empty toml should parse");
        assert_eq!(config.threads, None);
        assert_eq!(config.reference_files, None);
    }

    #[test]
    fn type_mismatches_caught_at_parse_time() {
        let path = temp_params_path("type-mismatch");
        fs::write(&path, r#"threads = "four""#).expect("write temp params file");

        let err = load_params_file(&path).expect_err("wrong field type should fail");
        assert!(err.to_string().contains("invalid type"));

        let _ = fs::remove_file(path);
    }
}
