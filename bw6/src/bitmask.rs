use ark_ff::{PrimeField, ToBytes, BigInteger};
use ark_std::convert::TryInto;

pub trait Bitmask<F: PrimeField>: ToBytes {
    fn to_bits(&self) -> Vec<bool>;

    fn to_chunks_as_field_elements(&self, bytes_in_chunk: usize) -> Vec<F>;
}

impl <F:PrimeField> Bitmask<F> for &[u8] {
    fn to_bits(&self) -> Vec<bool> {
        bytes_to_bits(&self)
    }

    fn to_chunks_as_field_elements(&self, bytes_in_chunk: usize) -> Vec<F> {
        self.chunks(bytes_in_chunk).map(bytes_to_field_elements::<F>).collect()
    }
}

fn bytes_to_field_elements<F: PrimeField>(chunk: &[u8]) -> F {
    let bytes_in_bigint =  F::BigInt::NUM_LIMBS * 8; // limbs are u64s
    let bytes_in_padding = bytes_in_bigint - chunk.len();
    let mut padding = vec![0u8; bytes_in_padding];
    let mut padded = Vec::with_capacity(bytes_in_bigint);
    padded.extend_from_slice(chunk);
    padded.append(&mut padding);
    let reversed: Vec<u8> = padded.iter().rev()
        .map(|byte|byte.reverse_bits())
        .collect();
    let limbs: Vec<u64> = reversed.rchunks_exact(8)
        .map(|limb| u64::from_be_bytes(limb.try_into().unwrap()))
        .collect();
    // let repr = F::BigInt::new(limbs.try_into().unwrap());
    let mut repr = F::BigInt::default();
    let repr_limbs = repr.as_mut();
    assert_eq!(repr_limbs.len(), limbs.len());
    repr_limbs.iter_mut().zip(limbs).for_each(|(a, b)| *a = b);
    F::from_repr(repr).unwrap()
}

fn byte_to_bits(byte: &u8) -> Vec<bool> {
    let mut res = vec![false; 8];
    let mut mask = 0b10000000u8;
    for i in 0..8 {
        if byte & mask == 1u8 {
            res[i] = true;
        }
        mask >>= 1;
    }
    res
}

fn bytes_to_bits(bytes: &[u8]) -> Vec<bool> {
    bytes.iter().map(|byte| byte_to_bits(byte)).flatten().collect()
}

#[cfg(test)]
mod tests {
    use crate::bitmask::bytes_to_field_elements;
    use ark_ff::{Zero, One, PrimeField, BigInteger};
    use ark_bw6_761::Fr;

    #[test]
    pub fn test_256_bit_chunks_for_bw6() {
        let bytes_in_chunk = 256 / 8;

        let mut bytes = vec![0u8; bytes_in_chunk]; // 0
        assert_eq!(bytes_to_field_elements::<Fr>(&bytes), Fr::zero());

        bytes[0] = 0b10000000u8; // 1
        assert_eq!(bytes_to_field_elements::<Fr>(&bytes), Fr::one());

        bytes[bytes_in_chunk - 1] = 0b00000001u8; // 2^255 + 1
        let one = <Fr as PrimeField>::BigInt::from(1u64);
        let mut bi = one.clone();
        bi.muln(255);
        bi.add_nocarry(&one);
        assert_eq!(bytes_to_field_elements::<Fr>(&bytes), Fr::from_repr(bi).unwrap());
    }
}