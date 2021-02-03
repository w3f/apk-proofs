use ark_bw6_761::{BW6_761, Fr};
use ark_ff::{FftField, FftParameters};
use crate::{kzg, KZG_BW6};
use ark_std::convert::TryInto;
use rand::Rng;


pub struct Setup {
    pub domain_size: usize,
    pub kzg_params: kzg::Params<BW6_761>,
}


impl Setup {
    pub fn generate<R: Rng>(log_domain_size: u32, rng: &mut R) -> Self {
        let domain_size = 2u64.pow(log_domain_size);
        let domain_size: usize = domain_size.try_into().expect("domain size doesn't fit usize");

        let max_poly_degree = 3 * domain_size - 3; // deg(q) = 3n-3
        let log_ceil_max_poly_degree = max_poly_degree.next_power_of_two().trailing_zeros();
        assert!(log_ceil_max_poly_degree <= <Fr as FftField>::FftParams::TWO_ADICITY, "");
        let kzg_params = KZG_BW6::setup(max_poly_degree, rng).unwrap();


        Self {
            domain_size,
            kzg_params,
        }
    }

    pub fn max_keyset_size(&self) -> u64 {
        (self.domain_size - 1) as u64 // the decrement accounts for the accumulator initial value
    }
}