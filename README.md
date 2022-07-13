# Fully Succinct BLS signature aggregation

Individual BLS signatures on the same message can be aggregated into a single signature
that can be verified in constant time, given the verifier knows the aggregate public key of the set of the actual signers [[1]](https://eprint.iacr.org/2018/483). 
However, computing the aggregate public key is linear in the number of the actual signers
and requires the verifier to know the individual public keys.

This repo contains PoC implementations of succinct proofs of correctness of the aggregate public key,
given the verifier knows the commitment to the list of public keys of all the eligible signers.

See [a code example](bw6/examples/recursive.rs) for a sketch of a blockchain light client design exploiting such proofs.

## Building the Specification document

The specification is built by means of [Cargo spec](https://crates.io/crates/cargo-spec) crate. To build the specification document, one can simply invoke:
```
$ cargo spec build
```
[`specification.md`](./specification.md) then shall contains the newly built specification. The specification could be easily converted to HTML by the help of `pandoc` if desired:

```
pandoc -f commonmark specification.md --standalone  --output specification.html 
```
