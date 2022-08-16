# Fully Succinct BLS signature aggregation

Individual BLS signatures on the same message can be aggregated into a single signature
that can be verified in constant time, given the verifier knows the aggregate public key of the set of the actual signers [[1]](https://eprint.iacr.org/2018/483). 
However, computing the aggregate public key is linear in the number of the actual signers
and requires the verifier to know the individual public keys.

This repo contains PoC implementations of succinct proofs of correctness of the aggregate public key,
given the verifier knows the commitment to the list of public keys of all the eligible signers.

See [a code example](bw6/examples/recursive.rs) for a sketch of a blockchain light client design exploiting such proofs.

[The formal description and security model for our succinct proofs as well as their application to accountable light clients can be found here](https://github.com/w3f/apk-proofs/blob/main/Light%20Client.pdf). [A video presentation of this work as part of ZK Summit 7 is available here](https://www.youtube.com/watch?v=UaPdDYarKGY&list=PLj80z0cJm8QFnY6VLVa84nr-21DNvjWH7&index=19).    
