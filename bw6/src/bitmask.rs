use ark_ff::{PrimeField, FpParameters, BitIteratorLE};
use ark_std::convert::{TryInto, TryFrom};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

const BITS_IN_LIMB: usize = 64;

/// A bitmask that can be encoded as a Vec of field elements.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Bitmask {
    limbs: Vec<u64>,
    padding_size: usize,
}

impl AsRef<[u64]> for Bitmask {
    fn as_ref(&self) -> &[u64] {
        &self.limbs
    }
}

// u64 limbs are used for internal representation for compatibility with ark_ff::BigInteger.
// Everything is little-endian, meaning that the least significant bit/limb is the leftmost.
// As a consequence 64-bit chunks are reversed to form equivalent u64 limbs.
// The highest limb of the bitmask is right padded with 0s.
impl Bitmask {

    pub fn from_bits(bits: &[bool]) -> Self {
        // repr = bitmask + padding
        let bitmask_size = bits.len();
        let limbs_required = div_ceil(bitmask_size, BITS_IN_LIMB);
        let repr_size = BITS_IN_LIMB * limbs_required;
        let padding_size = repr_size - bitmask_size;

        let mut repr_bits = Vec::with_capacity(repr_size);
        repr_bits.extend_from_slice(bits);
        repr_bits.extend_from_slice(&vec![false; padding_size]);

        let limbs_bits_iter = repr_bits.chunks_exact(BITS_IN_LIMB);
        assert_eq!(limbs_bits_iter.remainder().len(), 0);
        let limbs = limbs_bits_iter.map(|bits| bits_to_limb(bits.try_into().unwrap())).collect();
        Self { limbs, padding_size }
    }

    pub fn to_bits(&self) -> Vec<bool> {
        let mut bits: Vec<_> = BitIteratorLE::new(self).collect();
        bits.truncate(self.size()); // remove padding
        bits
    }

    pub fn to_bits_as_field_elements<F: PrimeField>(&self) -> Vec<F> {
        self.to_bits().iter()
            .map(|b| if *b { F::one() } else { F::zero() })
            .collect()
    }

    pub fn size(&self) -> usize {
        let repr_size = BITS_IN_LIMB * self.limbs.len();
        let bitmask_size = repr_size - self.padding_size;
        bitmask_size
    }

    // TODO: padding
    pub fn count_ones(&self) -> usize {
        self.limbs.iter().map(|limb| usize::try_from(limb.count_ones()).unwrap()).sum()
    }

    /// Splits the bits into chunks of the specified size and converts the chunks into field elements,
    /// interpreting the lowest bit of the chunk as the least significant (~little-endian).
    /// Panics if the chunk doesn't have the unique representation in the field (chunk size exceeds the field capacity).
    pub fn to_chunks_as_field_elements<F: PrimeField>(&self, limbs_in_chunk: usize) -> Vec<F> {
        let bits_in_chunk = BITS_IN_LIMB * limbs_in_chunk;
        assert!(bits_in_chunk <= F::Params::CAPACITY.try_into().unwrap());
        self.limbs.chunks(limbs_in_chunk).map(limbs_to_field_elements::<F>).collect()
    }
}


// may overflow
fn div_ceil(a: usize, b: usize) -> usize {
    (a + b - 1) / b
}

// Assembles a u64 limb from a little-endian bit slice.
fn bits_to_limb(bits: &[bool; 64]) -> u64 {
    bits.iter().rev().fold(0u64, |limb_acc, next_bit| (limb_acc << 1) ^ (*next_bit as u64))
}

fn limbs_to_field_elements<F: PrimeField>(limbs: &[u64]) -> F {
    let mut repr = F::BigInt::default();
    let repr_limbs = repr.as_mut();
    assert!(repr_limbs.len() >= limbs.len());
    repr_limbs.iter_mut().zip(limbs).for_each(|(a, b)| *a = *b);
    F::from_repr(repr).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bw6_761::Fr;
    use ark_std::test_rng;
    use crate::test_helpers::random_bits;

    pub fn _test_from_bits_to_bits(size: usize) {
        let bits = random_bits(size, 1.0 / 2.0, &mut test_rng());
        let bitmask = Bitmask::from_bits(&bits);
        let limbs_in_bitmask = div_ceil(size, BITS_IN_LIMB);
        assert_eq!(bitmask.limbs.len(), limbs_in_bitmask);
        assert_eq!(bitmask.padding_size, BITS_IN_LIMB * limbs_in_bitmask - size);
        assert_eq!(bitmask.size(), size);
        assert_eq!(bits, bitmask.to_bits())
    }

    #[test]
    pub fn test_bitmask_from_bits_to_bits() {
        _test_from_bits_to_bits(0);
        _test_from_bits_to_bits(63);
        _test_from_bits_to_bits(64);
        _test_from_bits_to_bits(65);
        _test_from_bits_to_bits(1024);
    }

    #[test]
    pub fn test_bitmask_limbs_to_field_elements() {
        let bits = random_bits(1024, 1.0 / 2.0, &mut test_rng());
        let bitmask = Bitmask::from_bits(&bits);
        let chunks = bitmask.to_chunks_as_field_elements::<Fr>(1);
        assert_eq!(chunks, bitmask.limbs.iter().cloned().map(Fr::from).collect::<Vec<_>>());
    }

    #[test]
    pub fn test_bitmask_bigger_chunks_to_field_elements() {
        let bits = random_bits(1024, 1.0 / 2.0, &mut test_rng());
        let bitmask = Bitmask::from_bits(&bits);
        let chunks = bitmask.to_chunks_as_field_elements::<Fr>(4);
        assert_eq!(chunks.len(), 4);
    }


    pub fn _test_to_field_element(set_bit_positions: Vec<usize>, expected: u128) {
        let mut bits = vec![false; 128];
        assert!(set_bit_positions.len() <= 128);
        for i in set_bit_positions {
            bits[i] = true;
        }

        let bitmask = Bitmask::from_bits(&bits);
        let chunks = bitmask.to_chunks_as_field_elements::<Fr>(2);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0], Fr::from(expected));
    }

    #[test]
    pub fn test_to_field_elements() {
        _test_to_field_element(vec![0], 1);
        _test_to_field_element(vec![63], 2u128.pow(63));
        _test_to_field_element(vec![0, 64], 2u128.pow(64) + 1);
        _test_to_field_element(vec![126], 2u128.pow(126));
    }
}