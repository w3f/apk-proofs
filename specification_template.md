# APK Proofs

APK proofs are SNARK based proof protocols to enable verification of aggregated BLS signatures without the knowledge of individual public keys of all the signers. APK proofs uses its own custom SNARK. The protocol is implemented in two flavours: [BW6](#BW6) and [MIPP](#MIPP) which are specified in the following sections.

## BW6
### Preliminaries
In this section we describe the preliminary data structure which are being used by APK prover and verifier. 

#### The Key Set
{sections.bw6-keyset}

### Prover

{sections.bw6-prover}
### BW6 Example
This section sketches a blockchain light client design exploiting APK proofs.

{sections.bw6-example}

## MIPP
