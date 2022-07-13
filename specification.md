# APK Proofs

APK proofs are SNARK based proof protocols to enable verification of aggregated BLS signatures without the knowledge of individual public keys of all the signers. APK proofs uses its own custom SNARK. The protocol is implemented in two flavours: [BW6](#BW6) and [MIPP](#MIPP) which are specified in the following sections.

## BW6
### Prover

Prover is responsible to generate APK proofs. The `Prover` struct encapsultes this tasks. It contains the following fields:
- `Domains`: ???
- `Keyset`: ???
- `KzgCommitterKey`: the set points in G1 corresponding to $\tau^n G_1$.
- `Transcript`: Representing the statement which purportedly has been signed by the aggregated public key.

### BW6 Example
This section sketches a blockchain light client design exploiting APK proofs.



## MIPP
