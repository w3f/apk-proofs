use ark_bw6_761::{BW6_761, Fr};
use ark_ff::FftField;
use fflonk::pcs::kzg::urs::URS;
use fflonk::pcs::PCS;
use rand::Rng;

use crate::NewKzgBw6;

pub fn generate_for_keyset<R: Rng>(keyset_size: usize, rng: &mut R) -> URS<BW6_761> {
    let required_domain_size = keyset_size + 1; // additional slot is occupied by affine addition accumulator initial value
    // as we use radix 2 domains
    let required_domain_size = required_domain_size.next_power_of_two();
    let log_domain_size = required_domain_size.trailing_zeros();
    generate_for_domain(log_domain_size, rng)
}

pub fn generate_for_domain<R: Rng>(log_domain_size: u32, rng: &mut R) -> URS<BW6_761> {
    let domain_size = 2usize.pow(log_domain_size);
    // to operate with polynomials of degree up to 4 * domain_size, there should exist a domain of size 4 * domain_size
    assert!(log_domain_size + 2 <= Fr::TWO_ADICITY, "not enough 2-adicity in ark_bw6_761::Fr");

    // the highest degree polynomial prover needs to commit is the quotient q=aggregate_constraint_polynomial/vanishing_polynomial
    // as the highest constraint degree is 4n-3, deg(q) = 3n-3
    let max_poly_degree = highest_degree_to_commit(domain_size);
    let kzg_params = NewKzgBw6::setup(max_poly_degree, rng);
    // assert!(kzg_params.fits(domain_size));
    kzg_params
}

fn highest_degree_to_commit(domain_size: usize) -> usize {
    3 * domain_size - 3
}

// impl kzg::Params<BW6_761> {
//     pub fn fits(&self, domain_size: usize) -> bool {
//         highest_degree_to_commit(domain_size) <= self.get_pk().max_degree()
//     }
// }