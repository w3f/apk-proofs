use std::env;

pub fn parse_args_or(def_log_n: u32, scheme: &str) -> u32 {
    let arg1 = env::args().skip(1).next();
    match arg1 {
        None => {
            println!("LOG_N parameter is not provided, using the default.\n\
            Run with '--test {} LOG_N', where LOG_N is the domain size binary log", scheme);
            def_log_n
        }
        Some(arg_log_n) => {
            arg_log_n.parse().expect(&format!("{} is not a valid parameter", arg_log_n))
        }
    }
}