# APK Proofs

APK proofs are SNARK based proof protocols to enable verification of aggregated BLS signatures without the knowledge of individual public keys of all the signers. APK proofs uses its own custom SNARK. The protocol is implemented in two flavours: [BW6](#BW6) and [MIPP](#MIPP) which are specified in the following sections.

## BW6
### Preliminaries
In this section we describe the preliminary data structure which are being used by APK prover and verifier. 

#### The Key Set
APK proof provides polynomial commitment to the vector of public keys. As such the fundamental structure used is the set of public //~ keys which prover has to commit to. The `Keyset` struct represent that set. Whereas `KeysetCommitment` struct is used to store
Let 'pks' be such a vector that commit(pks) == KeysetCommitment::pks_comm, also let
domain_size := KeysetCommitment::domain.size and
keyset_size := KeysetCommitment::keyset_size
Then the verifier needs to trust that:
1.a. pks.len() == KeysetCommitment::domain.size
  b. pks[i] lie in BLS12-377 G1 for i=0,...,domain_size-2
  c. for the 'real' keys pks[i], i=0,...,keyset_size-1, there exist proofs of possession
     for the padding, pks[i], i=keyset_size,...,domain_size-2, dlog is not known,
     e.g. pks[i] = hash_to_g1("something").
  pks[domain_size-1] is not a part of the relation (not constrained) and can be anything,
  we set pks[domain_size-1] = (0,0), not even a curve point.
2. KeysetCommitment::domain is the domain used to interpolate pks
In light client protocols the commitment is to the upcoming validator set, signed by the current validator set.
Honest validator checks the proofs of possession, interpolates with the right padding over the right domain,
computes the commitment using the right parameters, and then sign it.
Verifier checks the signatures and can trust that the properties hold under some "2/3 honest validators" assumption.
As every honest validator generates the same commitment, verifier needs to check only the aggregate signature.
Actual public keys, no padding.
Interpolations of the coordinate vectors of the public key vector WITH padding.
Domain used to compute the interpolations above.
Polynomials above, evaluated over a 4-times larger domain.
Used by the prover to populate the AIR execution trace.


### Prover

Prover is responsible to generate APK proofs. The `Prover` struct encapsultes this tasks. It contains the following fields:
- `Domains`: ???
- `Keyset`: ???
- `KzgCommitterKey`: the set points in G1 corresponding to $\tau^n G_1$.
- `Transcript`: Representing the statement which purportedly has been signed by the aggregated public key.

### BW6 Example
This section sketches a blockchain light client design exploiting APK proofs.



## MIPP
