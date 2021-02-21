# Fully Succinct BLS signature aggregation

Individual BLS signatures on the same message can be aggregated into a single signature
that can be verified in constant time, given the verifier knows the aggregate public key of the set of actual signers [[1]](https://eprint.iacr.org/2018/483). 
However, computing the aggregate public key is linear in the number of the actual signers
and requires the verifier to know the individual public keys.

This repo contains PoC implementations of succinct proofs of correctness of the aggregate public key,
given the verifier knows the commitment to the list of public keys of all the eligible signers.