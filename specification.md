# APK Proofs

APK proofs are SNARK based proof protocols to enable verification of aggregated BLS signatures without the knowledge of individual public keys of all the signers. APK proofs uses its own custom SNARK. The protocol is implemented in two flavours: [BW6](#BW6) and [MIPP](#MIPP) which are specified in the following sections.

## BW6
### Preliminaries
In this section we describe the preliminary data structure which are being used by APK prover and verifier. 

#### The Key Set
In short, APK proof provides polynomial commitment to the vector of public keys.
Furthermore, it offers the mean to use this commitment to verify BLS signatures signed
by the subset of those public keys. As a result, the verifier does not need to know the
public keys to verify aggregated BLS signature signed by them.

In light client protocols, such commitment is used to commit to the upcoming validator set, signed by the current validator set.
Honest validator should check the proofs of possession of each public key belong to an upcoming validator and arrange them in a
sequence with a determistic order. They then should deterministically pad the sequence to get a vector of consistent length
of the right domain which they use to interpolate the polynomial and to compute the commitment using the right parameters, and then sign it.

Verifier checks the signatures and can trust that the properties hold under some "2/3 honest validators" assumption.
As every honest validator generates the same commitment, verifier needs to check only the aggregate signature.

As such the fundamental structure used is the set of public
keys which prover has to commit to. The `Keyset` struct represent that set. Whereas `KeysetCommitment` struct is used to store
prover's commitment to the key set.

Let 'pks' be such a vector that

 commit(pks) == KeysetCommitment::pks_comm,

 also let:

`domain_size := KeysetCommitment::domain.size`

 and

`keyset_size := KeysetCommitment::keyset_size`

 Then the verifier needs to trust that:

1. The following:
    - `pks.len() == KeysetCommitment::domain.size`
    - `pks[i]` lie in BLS12-377 G1 for `i=0,...,domain_size-2`
    - For the 'real' `keys pks[i]`, `i=0,...,keyset_size-1`, there exist proofs of possession
    - For the padding, `pks[i], i=keyset_size,...,domain_size-2`, `dlog is not known,
      e.g. pks[i] = hash_to_g1("something").
     - `pks[domain_size-1]` is not a part of the relation (not constrained) and can be anything,
  we set pks[domain_size-1] = (0,0), not even a curve point.
2. `KeysetCommitment::domain` is the domain used to interpolate pks

*pks_comm*: Per-coordinate KZG commitments to a vector of BLS public keys on BLS12-377 represented in affine.
$$([pkx]({\tau}), [pky]({\tau})$$ where:

$$ pkx(X) = \sum_{i=0}^{n-1} pkx_i \cdot L_i(X). $$
$$ pky(X) = \sum_{i=0}^{n-1} pky_i \cdot L_i(X). $$
Domain used to interpolate the vectors above. Radix2 Domain Works only for fields
that have a large multiplicative subgroup of size that is a power-of-2.
The actual size of keyset i.e. the number of possible signers in contrast to the size of keyset vector after padding
Actual public keys in form of Projective Points on G1, no padding.
Interpolations of the coordinate vectors of the public key vector which includes dummy padded keys 
Domain used to compute the interpolations above.
Polynomials above, evaluated over a 4-times larger domain.
Used by the prover to populate the AIR execution trace.


### Prover

Prover is responsible to generate APK proofs. The `Prover` struct encapsultes this task. It contains the following fields:
- `Domains`: ???
- `Keyset`: set of all committe public keys (?)
- `KzgCommitterKey`: the set points in G1 corresponding to $\tau^n G_1$.
- `Transcript`: Representing the statement which purportedly has been signed by the aggregated public key.
Prover::new give the set of committees key, commitment to the set and the set of kzg points used for the
commitment, constructs a new prover. 

### BW6 Example
This section sketches a blockchain light client design exploiting APK proofs.



## MIPP