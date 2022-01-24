# Accountable BLS signature aggregation

## What is it all about

In a recent popular [post](https://moxie.org/2022/01/07/web3-first-impressions.html) Moxie Marlinspike, Signal creator, shares his concerns on the current state of web3. His major point is that web3 is still interacted with through centralized intermediates. W3F sees the solution to that in broader adoption of techniques allowing trustless communication within web3 without running a personal node -- light clients.

This repo provides utilities for building efficient cryptographic light clients for PoS blockchains with validators using BLS signatures to sign on consensus messages, such as validator set change or block finality (future Kusama/Polkadot).

## How it works

It introduces a concept of public 'chain' key that can be used to verify signatures collectively produced by the validators on behalf of the chain in efficient and accountable manner. The 'chain' key is not the BLS aggregate of the validators' public keys, usable only if the signature of each validator is available, that is impossible due to malicious signers and network issues. Neither it is a threshold public key, that's generation is hard to scale.

Instead, the 'chain' key is a succinct commitment to the set of active validators' BLS public keys. It allows a trustless helper aggregate the public keys of the actual signers and generate a proof that the aggregated keys are eligable. Moreover, the helper produces the bitmask, identifying who on the validator list signed. It also builds the aggregate signature and passes it to a light client together with the aggregate public key, the bitmask, and the proof of their validity.    

The light client uses the 'chain' key to verify the proof. After that, the key received from the helper can be used to check the signature. 

When the validator set changes, the new 'chain' key is signed by the current validators, so the current 'chain' key may be used to validate the new one in the same way. If a light client doesn't have the current 'chain' key it has to catch up validating the next 'chain' key against the previous one on every validator set change starting from the genesys. A newly joining light client might even have to start this recursive validation from the chain genesis. [A code example](bw6/examples/recursive.rs) has some additional explanations.  

## In other words
Effectively the novel construction implements fully succinct (avoids linear key aggregation), accountable (stronger than threshold) BLS signature aggregation. It doesn't require complex setup, while permitting keys rotation. As a core signature scheme common BLS with proofs of possessions is used.
