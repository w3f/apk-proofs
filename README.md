# Fully Succinct BLS Signature Aggregation

Individual BLS signatures on the same message can be aggregated into a single signature
that can be verified in constant time, given the verifier knows the aggregate public key of the set of the actual signers [[1]](https://eprint.iacr.org/2018/483). 
However, computing the aggregate public key is linear in the number of the actual signers and requires the verifier to know the individual public keys.

We avoid such heavy computation for verifiers that are constrained resource-wise and computation-wise (e.g., mobile phones, smart contracts on  blockchains) by desining custom non-interactive succinct arguments of knowledge (SNARKs) that compute and ensure the correctness of an apk, i.e., an aggregated public key of actual signers. This repo contains PoC implementations as well as formalisations for our custom SNARKs for apk, given the verifier knows only a commitment to the list of public keys of all the eligible signers and a bitmask identifying the actuall sigers of a message. 

See [a code example](bw6/examples/recursive.rs) for a sketch of a blockchain light client design exploiting such proofs.

# Formal Write-up
The formal description and security model for our custom succinct arguments as well as their application to accountable light clients for PoS blockchains can be found [here/stable version](https://eprint.iacr.org/2022/1205) and [here/on-going updates](https://github.com/w3f/apk-proofs/blob/main/Light%20Client.pdf). A high-level summary of this work as well as its connection with related reseach effort from the Web3 Foundation can be found [here](https://research.web3.foundation/en/latest/polkadot/LightClientsBridges/index.html).

# Video Presentations
A video presentation of this work at sub0 2022 is [available here](https://www.youtube.com/watch?v=MCvX9ZZhO4I&list=PLOyWqupZ-WGvywLqJDsMIYdCn8QEa2ShQ&index=19) ([slides](https://docs.google.com/presentation/d/16LlsXWY2Q6_6QGZxkg84evaJqWNk6szX)), and at ZK Summit 7 is [available here](https://www.youtube.com/watch?v=UaPdDYarKGY&list=PLj80z0cJm8QFnY6VLVa84nr-21DNvjWH7&index=19).    

# How to Reproduce Results in the Formal Write-up

1. [Install](https://www.rust-lang.org/tools/install) Rust toolchain.
2. Run one of the commands bellow.

#### Basic Accountable Scheme 
> cargo test --release --features "parallel print-trace" --test basic 10

> cargo test --release --features "parallel print-trace" --test basic 16 

> cargo test --release --features "parallel print-trace" --test basic 20
 
<br/>

#### Packed Accountable Scheme
> cargo test --release --features "parallel print-trace" --test packed 10

> cargo test --release --features "parallel print-trace" --test packed 16

> cargo test --release --features "parallel print-trace" --test packed 20

<br/>

#### Counting Scheme
> cargo test --release --features "parallel print-trace" --test counting 10

> cargo test --release --features "parallel print-trace" --test counting 16

> cargo test --release --features "parallel print-trace" --test counting 20

<br/>

The output should look, for example, like
```
Running test for the 'basic' scheme for N = 2^10
Start:   Prover
End:     Prover ....................................................................520.741ms
Start:   Verifier
End:     Verifier ..................................................................25.871ms
```
