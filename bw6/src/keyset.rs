use crate::{hash_to_curve, NewKzgBw6};
use ark_poly::{Radix2EvaluationDomain, EvaluationDomain, Evaluations};
use ark_bls12_377::G1Projective;
use ark_bw6_761::Fr;
use ark_poly::univariate::DensePolynomial;
use ark_ec::ProjectiveCurve;
use crate::domains::Domains;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};
use fflonk::pcs::kzg::params::KzgCommitterKey;
use fflonk::pcs::{PCS, CommitterKey};


//~ In short, APK proof provides polynomial commitment to the vector of public keys.
//~ Furthermore, it offers the mean to use this commitment to verify BLS signatures signed
//~ by the subset of those public keys. As a result, the verifier does not need to know the
//~ public keys to verify aggregated BLS signature signed by them.
//~
//~ In light client protocols, such commitment is used to commit to the upcoming validator set, signed by the current validator set.
//~ Honest validator should check the proofs of possession of each public key belong to an upcoming validator and arrange them in a
//~ sequence with a determistic order. They then should deterministically pad the sequence to get a vector of consistent length
//~ of the right domain which they use to interpolate the polynomial and to compute the commitment using the right parameters, and then sign it.
//~
//~ Verifier checks the signatures and can trust that the properties hold under some "2/3 honest validators" assumption.
//~ As every honest validator generates the same commitment, verifier needs to check only the aggregate signature.
//~
//~ As such the fundamental structure used is the set of public
//~ keys which prover has to commit to. The `Keyset` struct represent that set. Whereas `KeysetCommitment` struct is used to store
//~ prover's commitment to the key set.
//~


//~ Let 'pks' be such a vector that
//~
//~  commit(pks) == KeysetCommitment::pks_comm,
//~
//~  also let:
//~
//~ `domain_size := KeysetCommitment::domain.size`
//~
//~  and
//~
//~ `keyset_size := KeysetCommitment::keyset_size`
//~
//~  Then the verifier needs to trust that:
//~
//~ 1. The following:
//~     - `pks.len() == KeysetCommitment::domain.size`
//~     - `pks[i]` lie in BLS12-377 G1 for `i=0,...,domain_size-2`
//~     - For the 'real' `keys pks[i]`, `i=0,...,keyset_size-1`, there exist proofs of possession
//~     - For the padding, `pks[i], i=keyset_size,...,domain_size-2`, `dlog is not known,
//~       e.g. pks[i] = hash_to_g1("something").
//~      - `pks[domain_size-1]` is not a part of the relation (not constrained) and can be anything,
//~   we set pks[domain_size-1] = (0,0), not even a curve point.
//~ 2. `KeysetCommitment::domain` is the domain used to interpolate pks
//~
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct KeysetCommitment {
    //~ *pks_comm*: Per-coordinate KZG commitments to a vector of BLS public keys on BLS12-377 represented in affine.
    //~ $$([pkx]({\tau}), [pky]({\tau})$$ where:
    //~
    //~ $$pkx(X) = \sum_{i=0}^{n-1} pkx_i \cdot L_i(X).$$
    //~ $$pky(X) = \sum_{i=0}^{n-1} pky_i \cdot L_i(X).$$
    pub pks_comm: (ark_bw6_761::G1Affine, ark_bw6_761::G1Affine),
    //~ Domain used to interpolate the vectors above. Radix2 Domain Works only for fields
    //~ that have a large multiplicative subgroup of size that is a power-of-2.
    pub domain: Radix2EvaluationDomain<Fr>, // could be defined by it's generator
    //~ The actual size of keyset i.e. the number of possible signers in contrast to the size of keyset vector after padding
    pub keyset_size: usize
}

#[derive(Clone)]
pub struct Keyset {
    //~ Actual public keys in form of Projective Points on G1, no padding.
    pub pks: Vec<G1Projective>,
    //~ Interpolations of the coordinate vectors of the public key vector which includes dummy padded keys 
    pub pks_polys: [DensePolynomial<Fr>; 2],
    //~ Domain used to compute the interpolations above.
    pub domain: Radix2EvaluationDomain<Fr>,
    //~ Polynomials above, evaluated over a 4-times larger domain.
    //~ Used by the prover to populate the AIR execution trace.
    pub pks_evals_x4: Option<[Evaluations<Fr, Radix2EvaluationDomain<Fr>>; 2]>,
}

impl Keyset {
    pub fn new(pks: Vec<G1Projective>) -> Self {
        let min_domain_size = pks.len() + 1; // extra 1 accounts apk accumulator initial value
        let domain = Radix2EvaluationDomain::<Fr>::new(min_domain_size).unwrap();

        let mut padded_pks = pks.clone();
        // a point with unknown discrete log
        let padding_pk = hash_to_curve::<ark_bls12_377::G1Projective>(b"apk-proofs");
        padded_pks.resize(domain.size(), padding_pk);

        // convert into affine coordinates to commit
        let (pks_x, pks_y) = G1Projective::batch_normalization_into_affine(&padded_pks).iter()
            .map(|p| (p.x, p.y))
            .unzip();
        let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, domain).interpolate();
        let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, domain).interpolate();
        Self {
            pks,
            domain,
            pks_polys: [pks_x_poly, pks_y_poly],
            pks_evals_x4: None,
        }
    }

    // Actual number of signers, not including the padding
    pub fn size(&self) -> usize {
        self.pks.len()
    }

    pub fn amplify(&mut self) {
        let domains = Domains::new(self.domain.size());
        let pks_evals_x4 = self.pks_polys.clone().map(|z| domains.amplify_polynomial(&z));
        self.pks_evals_x4 = Some(pks_evals_x4);
    }

    pub fn commit(&self, kzg_pk: &KzgCommitterKey<ark_bw6_761::G1Affine>) -> KeysetCommitment {
        assert!(self.domain.size() <= kzg_pk.max_degree() + 1);
        let pks_x_comm= NewKzgBw6::commit(kzg_pk, &self.pks_polys[0]).0;
        let pks_y_comm= NewKzgBw6::commit(kzg_pk, &self.pks_polys[1]).0;
        KeysetCommitment {
            pks_comm: (pks_x_comm, pks_y_comm),
            domain: self.domain,
            keyset_size: self.pks.len(),
        }
    }

    pub fn aggregate(&self, bitmask: &[bool]) -> ark_bls12_377::G1Projective {
        assert_eq!(bitmask.len(), self.size());
        bitmask.iter()
            .zip(self.pks.iter())
            .filter(|(b, _p)| **b)
            .map(|(_b, p)| p)
            .sum()
    }
}
