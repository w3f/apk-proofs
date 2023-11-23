use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ec::pairing::Pairing;
use ark_ff::{FftField, One, Zero};
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain};
use ark_poly::univariate::DensePolynomial;
use ark_std::UniformRand;
use fflonk::pcs::kzg::urs::URS;
use fflonk::pcs::PcsParams;

struct SetupKey<E: Pairing> {
    // first N powers
    taus_g1_tail: Vec<E::G1Affine>,
    // last N powers
    taus_g1_head: Vec<E::G1Affine>,
    tau_m_g1: E::G1Affine,
    //TODO
    tau_m_g2: E::G2Affine,
    lis_g1: Vec<E::G1Affine>,
    lis_g2: Vec<E::G2Affine>,
    domain: Radix2EvaluationDomain<E::ScalarField>,
    tau_n_minus_1_g2: E::G2Affine,
    tau_g2: E::G2Affine,
    // TODO: prepared
    g1: E::G1Affine,
    g2: E::G2Affine,
}

struct OnChainPublicKey<E: Pairing> {
    // Participant's index.
    i: usize,
    // BLS public key in G1.
    pk_g1: E::G1Affine,
    // BLS public key in G2.
    pk_g2: E::G2Affine,
    // Commitment to p_i(X) = n.sk.L_i(X) in G1.
    p_g1: E::G1Affine,
    // Commitment to p_i(X) = n.sk.L_i(X) in G2.
    p_g2: E::G2Affine,
    // Commitment to r_i(X) = (p_i(X) - sk) / X in G1. This is the KZG proof of p_i(0) = sk.
    r_g1: E::G1Affine,
    // Commitment to X^{M-(n-1)}.p_i(X), where (tau^M).G1 is the highest power in KZG setup.
    s_g1: E::G1Affine,
}

struct Shares<E: Pairing> {
    /// Participant's index.
    i: usize,
    /// Commitments to W_ij for this i and all j's.
    wi: Vec<E::G1Affine>,
}

impl<E: Pairing> SetupKey<E> {
    fn from_srs(srs: URS<E>, log_n: u32) -> Self {
        let n = 1 << log_n;
        assert!(srs.powers_in_g2.len() >= n + 1);
        let domain = Radix2EvaluationDomain::new(n).unwrap();
        let g1 = srs.powers_in_g1[0];
        let g2 = srs.powers_in_g2[0];
        let tau_n_minus_1_g2 = (srs.powers_in_g2[n].into_group() - srs.powers_in_g2[0]).into_affine();

        let mut lis_g1 = srs.powers_in_g1.iter()
            .take(n)
            .map(|p| p.into_group())
            .collect();
        domain.ifft_in_place(&mut lis_g1);
        let lis_g1 = E::G1::normalize_batch(&lis_g1);

        let mut lis_g2 = srs.powers_in_g2.iter()
            .take(n)
            .map(|p| p.into_group())
            .collect();
        domain.ifft_in_place(&mut lis_g2);
        let lis_g2 = E::G2::normalize_batch(&lis_g2);
        Self {
            taus_g1_tail: srs.powers_in_g1[..n].to_vec(),
            taus_g1_head: srs.powers_in_g1[(srs.powers_in_g1.len() - n)..].to_vec(),
            tau_m_g1: srs.powers_in_g1[srs.powers_in_g1.len() - n],
            tau_m_g2: srs.powers_in_g2[srs.powers_in_g2.len() - n - 1],
            lis_g1,
            lis_g2,
            domain,
            tau_n_minus_1_g2,
            tau_g2: srs.powers_in_g2[1],
            g1,
            g2,
        }
    }

    fn gen_pk(&self, i: usize, sk: E::ScalarField) -> (OnChainPublicKey<E>, Shares<E>) {
        let pk_g1 = (self.g1 * sk).into_affine();
        let pk_g2 = (self.g2 * sk).into_affine();
        let n = self.domain.size_as_field_element;
        let n_sk = sk * n;
        let p_g1 = (self.lis_g1[i] * &n_sk).into_affine();
        let p_g2 = (self.lis_g2[i] * &n_sk).into_affine();
        let wi = self.shares(i, sk);
        let p = self.p(i, sk);
        let r_g1 = E::G1::msm(&self.taus_g1_tail[..self.domain.size() - 1], &p[1..]).unwrap().into_affine();
        let s_g1 = E::G1::msm(&self.taus_g1_head, &p).unwrap().into_affine();
        (OnChainPublicKey {
            i,
            pk_g1,
            pk_g2,
            p_g1,
            p_g2,
            r_g1,
            s_g1,
        }, Shares {
            i,
            wi,
        })
    }

    fn verify_pk(&self, pk: &OnChainPublicKey<E>) {
        assert_eq!(E::pairing(pk.pk_g1, self.g2), E::pairing(self.g1, pk.pk_g2), "BLS public keys don't match");
        assert_eq!(E::pairing(pk.p_g1, self.g2), E::pairing(self.g1, pk.p_g2), "Commitments to p_i don't match");
        assert_eq!(E::pairing(pk.p_g1, self.g2), E::pairing(self.lis_g1[pk.i] * self.domain.size_as_field_element, pk.pk_g2), "pk.p_g1 malformed");
        assert_eq!(E::pairing(pk.r_g1, self.tau_g2), E::pairing(pk.p_g1.into_group() - pk.pk_g1.into_group(), self.g2), "pk.r malformed");
        // assert_eq!(E::pairing(pk.s_g1, self.g2), E::pairing(self.tau_m_g1, pk.p_g2), "pk.s malformed");
        assert_eq!(E::pairing(pk.s_g1, self.g2), E::pairing(pk.p_g1, self.tau_m_g2), "pk.s malformed");
    }

    fn verify_shares(&self, pk: &OnChainPublicKey<E>, shares: &Shares<E>) {
        assert_eq!(pk.i, shares.i);
        for (j, wij) in shares.wi.iter().enumerate() {
            let lhs = E::pairing(wij, self.tau_n_minus_1_g2);
            let c = if j != pk.i {
                pk.p_g1
            } else {
                (pk.p_g1.into_group() - pk.pk_g1 * self.domain.size_as_field_element).into_affine()
            };
            let rhs = E::pairing(c, self.lis_g2[j]);
            assert_eq!(lhs, rhs)
        }
    }

    fn p(&self, i: usize, sk: E::ScalarField) -> DensePolynomial<E::ScalarField> {
        let n = self.domain.size();
        let mut p = vec![E::ScalarField::zero(); n];
        p[i] = self.domain.size_as_field_element * sk;
        self.domain.ifft_in_place(&mut p);
        DensePolynomial::from_coefficients_vec(p)
    }

    fn shares(&self, i: usize, sk: E::ScalarField) -> Vec<E::G1Affine> {
        let roots = roots(&self.domain);
        let wi = roots[i];
        let li = self.lis_g1[i];
        let mut denoms: Vec<E::ScalarField> = roots.iter()
            .enumerate()
            .map(|(j, wj)|
                if j != i {
                    wi - wj
                } else {
                    E::ScalarField::one() // isn't used
                })
            .collect();
        ark_ff::batch_inversion(&mut denoms);
        // compute Cqj for j != i
        let mut cqs: Vec<E::G1> = roots.into_iter()
            .zip(self.lis_g1.iter())
            .zip(denoms).
            enumerate()
            .map(|(j, ((wj, &lj), cj))| {
                if j != i {
                    let c = cj * sk;
                    li * (wj * &c) - lj * (wi * &c)
                } else {
                    E::G1::zero() // doesn't affect the sum
                }
            }).collect();
        let cqi: E::G1 = -cqs.iter().sum::<E::G1>();
        cqs[i] = cqi;
        let cqs = E::G1::normalize_batch(&cqs);
        cqs
    }

    fn verify(&self, apk_proof: &ApkProof<E>, ct: E::G1Affine) {
        let a = E::pairing(ct, apk_proof.b_g2);
        let b = E::pairing(apk_proof.cs, self.g2);
        let c = E::pairing(apk_proof.cw, self.tau_n_minus_1_g2);
        assert_eq!(a, b + c);
        assert_eq!(E::pairing(apk_proof.cs, self.tau_m_g2), E::pairing(apk_proof.deg, self.g2));
        let n = self.domain.size_as_field_element;
        assert_eq!(E::pairing(apk_proof.cs.into_group() - apk_proof.apk_g1, self.g2), E::pairing(apk_proof.q0, self.tau_g2));
    }
}

struct ApkProverKey<E: Pairing> {
    cqtjs: Vec<E::G1Affine>,
    pks: Vec<OnChainPublicKey<E>>,
    lis_g2: Vec<E::G2Affine>,
}

struct ApkProof<E: Pairing> {
    /// BLS aggregate public key in G1.
    apk_g1: E::G1Affine,
    /// Bitmask //TODO: all keys vs verified ones
    b_g2: E::G2Affine,
    /// Commitment to p_S(X) = sum[sk_i.L_i(X) for i in S].
    cs: E::G1Affine,
    cw: E::G1Affine,
    q0: E::G1Affine,
    deg: E::G1Affine,
}

impl<E: Pairing> ApkProverKey<E> {
    fn prove_apk(&self, signers: &[usize]) -> ApkProof<E> { //TODO: can have duplicates
        let apk_g1 = signers.iter().map(|&i| self.pks[i].pk_g1).sum::<E::G1>().into();
        let b_g2 = signers.iter().map(|&i| self.lis_g2[i]).sum::<E::G2>().into();
        let cs = signers.iter().map(|&i| self.pks[i].p_g1).sum::<E::G1>().into();
        let q0 = signers.iter().map(|&i| self.pks[i].r_g1).sum::<E::G1>().into();
        let cw = signers.iter().map(|&i| self.cqtjs[i]).sum::<E::G1>().into();
        let deg = signers.iter().map(|&i| self.pks[i].s_g1).sum::<E::G1>().into();
        ApkProof {
            apk_g1,
            b_g2,
            cs,
            cw,
            q0,
            deg,
        }
    }
}

fn roots<F: FftField, D: EvaluationDomain<F>>(domain: &D) -> Vec<F> {
    let n = domain.size();
    let mut res = Vec::with_capacity(n);
    let mut wi = F::one();
    res.push(wi);
    for _ in 1..n {
        wi *= domain.group_gen();
        res.push(wi);
    }
    res
}

fn aggregate_shares<E: Pairing>(shares: &[Shares<E>]) -> Vec<E::G1Affine> {
    let n = shares.len();
    let mut qis = vec![vec![E::G1Affine::default(); n]; n];
    for j in 0..n {
        qis[j] = shares.iter()
            .map(|share_i| share_i.wi[j])
            .collect();
    }
    qis.iter()
        .map(|qi| qi.iter().sum::<E::G1>().into())
        .collect()
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective};
    use ark_std::{test_rng, UniformRand};

    use super::*;

    #[test]
    fn it_works() {
        let rng = &mut test_rng();

        let log_n = 2;
        let n = 1 << log_n;
        let srs = URS::<Bls12_381>::generate(n, n + 1, rng);
        let setup_key = SetupKey::from_srs(srs, log_n);

        let (pks, shares): (Vec<_>, Vec<_>) = (0..n)
            .map(|i| setup_key.gen_pk(i, Fr::rand(rng)))
            .unzip();

        for (pk_i, shares_i) in pks.iter()
            .zip(shares.iter()) {
            setup_key.verify_pk(pk_i);
            setup_key.verify_shares(pk_i, shares_i);
        }

        let agg_shares = aggregate_shares(&shares);

        let ct: G1Affine = pks.iter().map(|pk| pk.p_g1).sum::<G1Projective>().into();
        let prover_key = ApkProverKey {
            cqtjs: agg_shares,
            pks,
            lis_g2: setup_key.lis_g2.clone(),
        };

        let apk_proof = prover_key.prove_apk(&[1, 2]);
        setup_key.verify(&apk_proof, ct);

        let apk_proof = prover_key.prove_apk(&[3]);
        setup_key.verify(&apk_proof, ct);
    }
}
