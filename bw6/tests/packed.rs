mod helper;

fn main() {
    let log_n = helper::parse_args_or(8, "packed");
    println!("Running test for the 'packed' scheme for N = 2^{}", log_n);
    apk_proofs::test_helpers::test_packed_scheme(log_n);
}