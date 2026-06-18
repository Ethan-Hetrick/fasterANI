fasterANI: A reimplementation of FastANI in Rust, with better performance and optimized search modes.

```bash
# Build
export RUSTFLAGS="-C target-cpu=native"
cargo build --release

# Help menu
./target/release/fasterANI --help

# Run test data
./target/release/fasterANI 
    --reference assets/test-data/Escherichia_coli_str_K12_MG1655.fna \
    --query assets/test-data/Shigella_flexneri_2a_01.fna 2> /dev/null
```
