mod helper;

fn main() {
    let log_n = helper::parse_args_or(8, "counting");
    println!("Running test for the 'counting' scheme for N = 2^{}", log_n);
    apk_proofs::test_helpers::test_counting_scheme(log_n);
}