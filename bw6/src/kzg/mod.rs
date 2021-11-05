//! Here we construct a polynomial commitment that enables users to commit to a
//! single polynomial `p`, and then later provide an evaluation proof that
//! convinces verifiers that a claimed value `v` is the true evaluation of `p`
//! at a chosen point `x`. Our construction follows the template of the construction
//! proposed by Kate, Zaverucha, and Goldberg ([KZG11](http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf)).
//! This construction achieves extractability in the algebraic group model (AGM).

use ark_ec::msm::{FixedBaseMSM, VariableBaseMSM};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{One, PrimeField, UniformRand, Zero, FftField, Field};
use ark_poly::UVPolynomial;
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};

use ark_std::{format, marker::PhantomData, ops::Div, vec};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

use ark_std::{end_timer, start_timer};

use crate::utils;

// #[cfg(feature = "parallel")]
// use rayon::prelude::*;

use rand::RngCore;
use ark_std::ops::Mul;
use std::collections::HashSet;
use merlin::Transcript;
use crate::transcript::ApkTranscript;
use std::ops::Sub;

pub trait ShplonkTranscript<E: PairingEngine> {
    fn get_gamma(&mut self) -> E::Fr;
    fn commit_to_q(&mut self, point: &E::G1Affine);
    fn get_zeta(&mut self) -> E::Fr;
}

impl<E: PairingEngine> ShplonkTranscript<E> for Transcript {
    fn get_gamma(&mut self) -> E::Fr {
        let mut buf = [0u8; 16];
        self.challenge_bytes(b"gamma", &mut buf);
        E::Fr::from_random_bytes(&buf).unwrap()
    }

    fn commit_to_q(&mut self, point: &E::G1Affine) {
        let mut buf = vec![0; point.serialized_size()];
        point.serialize(&mut buf);
        self.append_message(b"q", &buf);
    }

    fn get_zeta(&mut self) -> E::Fr {
        let mut buf = [0u8; 16]; //TODO: full width
        self.challenge_bytes(b"zeta", &mut buf);
        E::Fr::from_random_bytes(&buf).unwrap()
    }
}

/// `Params` are the parameters for the KZG10 scheme.
#[derive(Clone, Debug)]
pub struct Params<E: PairingEngine> {
    /// Group elements of the form `{ \beta^i G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_g: Vec<E::G1Affine>,
    // TODO
    pub powers_of_g2: Vec<E::G2Affine>,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
}

impl<E: PairingEngine> Params<E> {
    pub fn get_pk(&self) -> ProverKey<E> {
        ProverKey(self.powers_of_g.clone()) //TODO: avoid cloning
    }

    pub fn get_vk(&self) -> VerifierKey<E> {
        VerifierKey {
            g: self.powers_of_g[0],
            h: self.h,
            beta_h: self.beta_h,
        }
    }

    pub fn get_evk(&self) -> ExtendedVerifierKey<E> {
        ExtendedVerifierKey {
            powers_in_g1: self.powers_of_g.clone(),
            powers_in_g2: self.powers_of_g2.clone(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct ProverKey<E: PairingEngine> (
    /// Group elements of the form `{ \beta^i G }`, where `i` ranges from 0 to `degree`.
    Vec<E::G1Affine>
);

impl<E: PairingEngine> ProverKey<E> {
    pub fn max_coeffs(&self) -> usize {
        self.0.len()
    }

    pub fn max_degree(&self) -> usize {
        self.max_coeffs() - 1
    }
}

#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerifierKey<E: PairingEngine> {
    /// The generator of G1.
    pub g: E::G1Affine,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
}

impl<E: PairingEngine> VerifierKey<E> {
    pub fn prepare(&self) -> PreparedVerifierKey<E> {
        PreparedVerifierKey {
            g: self.g,
            prepared_h: self.h.into(),
            prepared_beta_h: self.beta_h.into(),
        }
    }
}

/// `PreparedVerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(Clone, Debug)]
pub struct PreparedVerifierKey<E: PairingEngine> {
    /// The generator of G1.
    pub g: E::G1Affine,
    /// The generator of G2, prepared for use in pairings.
    pub prepared_h: E::G2Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    pub prepared_beta_h: E::G2Prepared,
}

#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
// Key for verifying n-point KZG with the polynomial and the proof committed to G1.
// Contains n powers in G1 and n+1 power in G2.
pub struct ExtendedVerifierKey<E: PairingEngine> {
    pub powers_in_g1: Vec<E::G1Affine>,
    pub powers_in_g2: Vec<E::G2Affine>,
}

/// `KZG10` is an implementation of the polynomial commitment scheme of
/// [Kate, Zaverucha and Goldbgerg][kzg10]
///
/// [kzg10]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
pub struct KZG10<E: PairingEngine, P: UVPolynomial<E::Fr>> {
    _engine: PhantomData<E>,
    _poly: PhantomData<P>,
}

impl<E, P> KZG10<E, P>
    where
        E: PairingEngine,
        P: UVPolynomial<E::Fr, Point=E::Fr> + From<DensePolynomial<E::Fr>>,
        for<'a, 'b> &'a P: Div<&'b P, Output=P>,
        for<'a, 'b> &'a P: Mul<&'b P, Output=P>,
        for<'a, 'b> &'a P: Sub<&'b P, Output=P>,
        for<'a> &'a P: Mul<E::Fr, Output=P>,
        for<'a> &'a P: Into<DenseOrSparsePolynomial<'a, E::Fr>>,

{
    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    pub fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Params<E> {
        let setup_time = start_timer!(|| format!("KZG10::Setup with degree {}", max_degree));
        let beta = E::Fr::rand(rng);
        let g = E::G1Projective::rand(rng);
        let h = E::G2Projective::rand(rng);

        let powers_of_beta = utils::powers(beta, max_degree);

        let window_size = FixedBaseMSM::get_mul_window_size(max_degree + 1);

        let scalar_bits = E::Fr::size_in_bits();
        let g_time = start_timer!(|| "Generating powers of G");
        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g);
        let powers_of_g = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
            scalar_bits,
            window_size,
            &g_table,
            &powers_of_beta,
        );
        end_timer!(g_time);

        let powers_of_g = E::G1Projective::batch_normalization_into_affine(&powers_of_g);


        let max_degree = 10; //TODO
        let powers_of_beta = utils::powers(beta, max_degree);
        let window_size = FixedBaseMSM::get_mul_window_size(max_degree + 1);
        let scalar_bits = E::Fr::size_in_bits();
        let g_time = start_timer!(|| "Generating powers of G");
        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, h);
        let powers_of_g2 = FixedBaseMSM::multi_scalar_mul::<E::G2Projective>(
            scalar_bits,
            window_size,
            &g_table,
            &powers_of_beta,
        );
        end_timer!(g_time);
        let powers_of_g2 = E::G2Projective::batch_normalization_into_affine(&powers_of_g2);


        let h = h.into_affine();
        let beta_h = h.mul(beta).into_affine();

        let pp = Params {
            powers_of_g,
            powers_of_g2,
            h,
            beta_h,
        };
        end_timer!(setup_time);
        pp
    }

    /// Outputs a commitment to `polynomial`.
    pub fn commit(
        powers: &ProverKey<E>,
        polynomial: &P,
    ) -> E::G1Affine {
        assert!(polynomial.degree() <= powers.max_degree());

        let commit_time = start_timer!(|| format!(
            "Committing to polynomial of degree {}",
            polynomial.degree()
        ));

        let (num_leading_zeros, plain_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(polynomial);

        let msm_time = start_timer!(|| "MSM to compute commitment to plaintext poly");
        let commitment = VariableBaseMSM::multi_scalar_mul(
            &powers.0[num_leading_zeros..],
            &plain_coeffs,
        );
        end_timer!(msm_time);
        end_timer!(commit_time);
        commitment.into()
    }

    /// Outputs a commitment to `polynomial`.
    pub fn commit_g<G: AffineCurve, PG: UVPolynomial<G::ScalarField>>(
        powers: &[G],
        polynomial: &PG,
    ) -> G
    {
        assert!(polynomial.degree() + 1 <= powers.len());

        let (num_leading_zeros, plain_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(polynomial);

        let commitment = VariableBaseMSM::multi_scalar_mul(
            &powers[num_leading_zeros..],
            &plain_coeffs,
        );

        commitment.into()
    }

    /// Compute witness polynomial.
    ///
    /// The witness polynomial w(x) the quotient of the division (p(x) - p(z)) / (x - z)
    /// Observe that this quotient does not change with z because
    /// p(z) is the remainder term. We can therefore omit p(z) when computing the quotient.
    pub fn compute_witness_polynomial(
        p: &P,
        point: P::Point,
    ) -> P {
        let divisor = P::from_coefficients_vec(vec![-point, E::Fr::one()]);

        let witness_time = start_timer!(|| "Computing witness polynomial");
        let witness_polynomial = p / &divisor;
        end_timer!(witness_time);

        witness_polynomial
    }

    pub(crate) fn open_with_witness_polynomial(
        powers: &ProverKey<E>,
        witness_polynomial: &P,
    ) -> E::G1Affine {
        assert!(witness_polynomial.degree() <= powers.max_degree());
        let (num_leading_zeros, witness_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(witness_polynomial);

        let witness_comm_time = start_timer!(|| "Computing commitment to witness polynomial");
        let w = VariableBaseMSM::multi_scalar_mul(
            &powers.0[num_leading_zeros..],
            &witness_coeffs,
        );
        end_timer!(witness_comm_time);

        w.into_affine()
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    pub fn open<'a>(
        powers: &ProverKey<E>,
        p: &P,
        point: P::Point,
    ) -> E::G1Affine {
        assert!(p.degree() <= powers.max_degree());
        let open_time = start_timer!(|| format!("Opening polynomial of degree {}", p.degree()));

        let witness_time = start_timer!(|| "Computing witness polynomials");
        let witness_poly = Self::compute_witness_polynomial(p, point);
        end_timer!(witness_time);

        let proof = Self::open_with_witness_polynomial(
            powers,
            &witness_poly,
        );

        end_timer!(open_time);
        proof
    }

    /// Verifies that `value` is the evaluation at `point` of the polynomial
    /// committed inside `comm`.
    pub fn check(
        pvk: &PreparedVerifierKey<E>,
        comm: &E::G1Affine,
        point: E::Fr,
        value: E::Fr,
        proof: E::G1Affine,
    ) -> bool {
        let t_check = start_timer!(|| "1-point KZG verification");

        let lhs = comm.into_projective() - &pvk.g.mul(value) + &proof.mul(point); // $[p(x)]_1 - y_1[1]_1 + x_1 [q(x)]_1$
        let pp = E::product_of_pairings(&[
            (lhs.into_affine().into(), pvk.prepared_h.clone()), // $e([p(x)]_1 - y_1[1]_1 + x_1 [q(x)]_1, [1]_2)$
            ((-proof).into(), pvk.prepared_beta_h.clone()), // $e([q(x)]_1,[x]_2)$
        ]);
        let valid = pp.is_one();
        end_timer!(t_check, || format!("valid = {}", valid));
        return valid;
    }

    pub fn aggregate_openings<R: RngCore>(
        vk: &PreparedVerifierKey<E>,
        commitments: &[E::G1Affine],
        points: &[E::Fr],
        values: &[E::Fr],
        proofs: &[E::G1Affine],
        rng: &mut R,
    ) -> (E::G1Projective, E::G1Projective) {
        let mut total_c = <E::G1Projective>::zero();
        let mut total_w = <E::G1Projective>::zero();

        let combination_time = start_timer!(|| "Combining commitments and proofs");
        let mut randomizer = E::Fr::one();
        // Instead of multiplying g in each turn, we simply accumulate
        // it's coefficients and perform a final multiplication at the end.
        let mut g_multiplier = E::Fr::zero();
        for (((c, z), v), w) in commitments.iter()
            .zip(points)
            .zip(values)
            .zip(proofs) {
            let mut temp = w.mul(*z); // $x_i [q_i(x)]_1$
            temp.add_assign_mixed(&c); // $[p_i(x)]_1 + x_i [q_i(x)]_1$
            let c = temp;
            g_multiplier += &(randomizer * v); // $r_i y_i$
            total_c += &c.mul(randomizer.into_repr()); // $r_i [p_i(x)]_1 + r_i x_i [q_i(x)]_1$
            total_w += &w.mul(randomizer); //  $r_i [q_i(x)]_1$
            // We don't need to sample randomizers from the full field,
            // only from 128-bit strings.
            randomizer = u128::rand(rng).into();
        }
        total_c -= &vk.g.mul(g_multiplier);
        // $(\sum_i r_i y_i) [1]_1$
        end_timer!(combination_time);

        (total_c, total_w)
    }

    pub fn batch_check_aggregated(
        vk: &PreparedVerifierKey<E>,
        total_c: E::G1Projective,
        total_w: E::G1Projective,
    ) -> bool {
        let to_affine_time = start_timer!(|| "Converting results to affine for pairing");
        let affine_points = E::G1Projective::batch_normalization_into_affine(&[total_c, -total_w]);
        let (total_c, total_w) = (affine_points[0], affine_points[1]);
        end_timer!(to_affine_time);

        let pairing_time = start_timer!(|| "Performing product of pairings");
        let valid = E::product_of_pairings(&[
            (total_c.into(), vk.prepared_h.clone()),
            (total_w.into(), vk.prepared_beta_h.clone()),
        ])
            .is_one();
        end_timer!(pairing_time);
        valid
    }

    /// Check that each `proof_i` in `proofs` is a valid proof of evaluation for
    /// `commitment_i` at `point_i`.
    pub fn batch_check<R: RngCore>(
        vk: &PreparedVerifierKey<E>,
        commitments: &[E::G1Affine],
        points: &[E::Fr],
        values: &[E::Fr],
        proofs: &[E::G1Affine],
        rng: &mut R,
    ) -> bool {
        let check_time =
            start_timer!(|| format!("Checking {} evaluation proofs", commitments.len()));
        let (total_c, total_w) = Self::aggregate_openings(vk, commitments, points, values, proofs, rng);
        let valid = Self::batch_check_aggregated(vk, total_c, total_w);
        end_timer!(check_time, || format!("Result: {}", valid));
        valid
    }

    /// Computes p = p0 + (r * p1) + (r^2 * p2) + ... + (r^n * pn)
    pub fn aggregate_polynomials(
        randomizer: E::Fr,
        polynomials: &[P],
    ) -> P {
        utils::randomize(randomizer, polynomials)
    }

    /// Computes C = C0 + (r * C1) + (r^2 * C2) + ... + (r^n * Cn)
    pub fn aggregate_commitments(
        randomizer: E::Fr,
        commitments: &[E::G1Affine],
    ) -> E::G1Affine {
        utils::horner(commitments, randomizer)
    }

    /// Computes v = v0 + (r * v1) + (r^2 * v2) + ... + (r^n * vn)
    pub fn aggregate_values(
        randomizer: E::Fr,
        values: &[E::Fr],
    ) -> E::Fr {
        utils::horner_field(values, randomizer)
    }

    // vanishing polynomial of the set
    // TODO: product set
    pub fn z_of_set<'a>(xs: impl IntoIterator<Item=&'a P::Point>) -> P {
        xs.into_iter()
            .map(|&xi| Self::z_of_point(xi))
            .reduce(|p, pi| &p * &pi)
            .unwrap()
    }

    pub fn shplonk_open_1_poly_n_points_interactive<T: ShplonkTranscript<E>>(
        powers: &ProverKey<E>,
        f: &P,
        xs: &HashSet<P::Point>,
        transcript: &mut T,
    ) -> (E::G1Affine, E::G1Affine) {
        assert!(f.degree() <= powers.max_degree());
        let z = Self::z_of_set(xs);
        let (q, r) = Self::divide_with_remainder(f, &z);
        let q1 = Self::open_with_witness_polynomial(powers, &q);

        transcript.commit_to_q(&q1);
        let zeta = transcript.get_zeta();

        let r_at_zeta = r.evaluate(&zeta);
        let z_at_zeta = z.evaluate(&zeta);

        let mut l = P::zero();
        l += f;
        l -= &(&q * z_at_zeta);
        l -= &P::from_coefficients_vec(vec![r_at_zeta]);

        let z_of_zeta = Self::z_of_point(zeta);
        let (q_of_l, r_of_l) = Self::divide_with_remainder(&l, &z_of_zeta);
        assert!(r_of_l.is_zero());
        let q_of_l1 = Self::open_with_witness_polynomial(powers, &q_of_l);
        (q1, q_of_l1)
    }

    pub fn shplonk_verify_1_poly_n_points_interactive<T: ShplonkTranscript<E>>(
        pvk: &PreparedVerifierKey<E>,
        f1: &E::G1Affine,
        q1: E::G1Affine,
        q_of_l1: E::G1Affine,
        xs: &[E::Fr],
        ys: &[E::Fr],
        transcript: &mut T,
    ) -> bool {
        transcript.commit_to_q(&q1);
        let zeta = transcript.get_zeta();

        let r = Self::interpolate(xs, ys);
        let r_zeta = r.evaluate(&zeta);
        let r_zeta_1 = pvk.g.mul(r_zeta);

        let z = Self::z_of_set(xs);
        let z_at_zeta = z.evaluate(&zeta);

        let l_1 = f1.into_projective() - r_zeta_1 - q1.mul(z_at_zeta) + q_of_l1.mul(zeta);

        E::product_of_pairings(&[
            (l_1.into_affine().into(), pvk.prepared_h.clone()),
            ((-q_of_l1).into(), pvk.prepared_beta_h.clone()),
        ]).is_one()
    }

    pub fn shplonk_open<T: ShplonkTranscript<E>>(
        powers: &ProverKey<E>,
        fs: &[P],
        xss: &[HashSet<P::Point>],
        transcript: &mut T,
    ) -> (E::G1Affine, E::G1Affine) {
        assert_eq!(xss.len(), fs.len(), "{} opening sets specified for {} polynomials", xss.len(), fs.len());
        let mut opening_set = HashSet::new();
        for xs in xss {
            opening_set.extend(xs);
        }
        let z = Self::z_of_set(&opening_set);

        let zs: Vec<P> = xss.iter().map(|xs| Self::z_of_set(xs)).collect();

        let (qs, rs): (Vec<_>, Vec<_>) = fs.iter().zip(&zs).map(|(fi, zi)| Self::divide_with_remainder(fi, zi)).unzip();

        let gamma = transcript.get_gamma();

        let q = utils::randomize(gamma, &qs);
        let q1 = Self::open_with_witness_polynomial(powers, &q);

        transcript.commit_to_q(&q1);
        let zeta = transcript.get_zeta();

        let z_zeta = z.evaluate(&zeta);
        let mut zs_zeta: Vec<_> = zs.iter().map(|zi| zi.evaluate(&zeta)).collect();
        let rs_zeta: Vec<_> = rs.iter().map(|ri| ri.evaluate(&zeta)).collect();
        ark_ff::batch_inversion(&mut zs_zeta);

        let gs = utils::powers(gamma, fs.len() - 1);

        let mut l = P::zero();
        for (((fi, ri), zi_inv), gi) in fs.iter()
            .zip(rs_zeta)
            .zip(zs_zeta)
            .zip(gs) {
            l += (gi * zi_inv, &(fi - &P::from_coefficients_vec(vec![ri])));
        }
        let l = &(&l - &q) * z_zeta;

        let z_of_zeta = Self::z_of_point(zeta);
        let (q_of_l, r_of_l) = Self::divide_with_remainder(&l, &z_of_zeta);
        assert!(r_of_l.is_zero());
        let q_of_l1 = Self::open_with_witness_polynomial(powers, &q_of_l);
        (q1, q_of_l1)
    }

    pub fn shplonk_verify<T: ShplonkTranscript<E>>(
        pvk: &PreparedVerifierKey<E>,
        f1s: &[E::G1Affine],
        q1: E::G1Affine,
        q_of_l1: E::G1Affine,
        xss: &Vec<Vec<E::Fr>>,
        yss: &Vec<Vec<E::Fr>>,
        transcript: &mut T,
    ) -> bool {
        let gamma = transcript.get_gamma();
        transcript.commit_to_q(&q1);
        let zeta = transcript.get_zeta();

        let mut opening_set = HashSet::new();
        for xs in xss {
            opening_set.extend(xs);
        }

        let z_at_zeta = opening_set.iter()
            .map(|xi| zeta - xi)
            .reduce(|z, zi| z * zi)
            .unwrap();


        let rs = xss.iter().zip(yss).map(|(xs, ys)| Self::interpolate(&xs, &ys));
        let rs_at_zeta = rs.map(|ri| ri.evaluate(&zeta));

        let mut zs_at_zeta: Vec<_> = xss.iter().map(|xs|
            xs.iter()
                .map(|xi| zeta - xi)
                .reduce(|z, zi| z * zi)
                .unwrap()
        ).collect();

        ark_ff::batch_inversion(&mut zs_at_zeta);

        let gs = utils::powers(gamma, f1s.len() - 1);

        let gzs: Vec<_> = gs.iter().zip(zs_at_zeta).map(|(&gi, zi_inv)| gi * zi_inv).collect();

        let f_1: E::G1Projective = f1s.iter().zip(&gzs).map(|(f1, gzi)| f1.mul(z_at_zeta * gzi)).sum();

        let r = rs_at_zeta.zip(gzs).map(|(ri_at_zeta, gzi)| ri_at_zeta * &gzi).sum::<E::Fr>() * z_at_zeta;

        let l_1 = f_1 - pvk.g.mul(r) - q1.mul(z_at_zeta) + q_of_l1.mul(zeta);

        E::product_of_pairings(&[
            (l_1.into_affine().into(), pvk.prepared_h.clone()),
            ((-q_of_l1).into(), pvk.prepared_beta_h.clone()),
        ]).is_one()
    }

    /// Flattens the matrix `fs` given as a list of rows by concatenating it's columns.
    /// Rows are right-padded by 0s to `max_row_len`.
    fn flatten<'a>(fs: impl ExactSizeIterator<Item=&'a [<E as PairingEngine>::Fr]>, max_row_len: usize) -> Vec<E::Fr> {
        let t = fs.len();
        let mut g = vec![E::Fr::zero(); t * max_row_len];
        for (i, fi) in fs.enumerate() {
            for (j, fij) in fi.iter().enumerate() {
                g[t * j + i] = *fij;
            }
        }
        g
    }

    /// `g = combine(f1,..,ft)`
    fn combine(fs: &[P]) -> P {
        let max_coeffs_len = fs.iter().map(|fi| fi.coeffs().len()).max().expect("empty slice");
        let fs = fs.iter().map(|fi| fi.coeffs());
        let g = Self::flatten(fs, max_coeffs_len);
        P::from_coefficients_vec(g)
    }

    /// On input a polynomial `g = combine(f1,...,ft)` and a point `x`,
    /// outputs single proof for `f1(x),...,ft(x)`.
    pub fn open_combined(
        powers: &ProverKey<E>,
        g: &P,
        x: P::Point,
        t: usize,
    ) -> E::G1Affine {
        assert!(g.degree() <= powers.max_degree()); //TODO:
        let q = Self::fflonk_quotient_poly(g, t, x);
        Self::open_with_witness_polynomial(powers, &q)
    }

    pub fn flonk_n_point_z(xs: &[P::Point], t: usize) -> P {
        xs.iter()
            .map(|&xi| Self::fflonk_vanishing_poly(t, xi))
            .reduce(|p, pi| &p * &pi)
            .unwrap()
    }

    pub fn flonk_open_1_poly_n_points(
        powers: &ProverKey<E>,
        g: &P,
        xs: &[P::Point],
        t: usize,
    ) -> E::G1Affine {
        assert!(g.degree() <= powers.max_degree()); //TODO:
        let z = Self::flonk_n_point_z(xs, t);
        let q = g / &z; // ignores the remainder
        Self::open_with_witness_polynomial(powers, &q)
    }

    /// (fflonk) opening polynomial for `g = combine(f1,...,ft)` in point `x`
    /// `q = g // X^t - x`
    fn fflonk_quotient_poly(
        g: &P,
        t: usize,
        x: P::Point,
    ) -> P {
        g / &Self::fflonk_vanishing_poly(t, x) // ignores the remainder
    }

    /// Primitive root of unity of order `t`.
    fn primitive_root(t: usize) -> E::Fr {
        E::Fr::get_root_of_unity(t).expect("no root of unity")
    }

    /// `z, zw,...,zw^{t-1}`
    fn fflonk_roots(z: E::Fr, t: usize) -> Vec<E::Fr> {
        let w = Self::primitive_root(t);
        let mut acc = z;
        let mut roots = vec![z];
        roots.resize_with(t, || {
            acc *= w;
            acc
        });
        roots
    }

    // Z(X) = X - x
    fn z_of_point(x: P::Point) -> P {
        P::from_coefficients_vec(vec![-x, E::Fr::one()])
    }

    // Z(X) = X^t - x
    fn fflonk_vanishing_poly(t: usize, x: P::Point) -> P {
        let mut z = vec![E::Fr::zero(); t + 1]; // deg(Z) = t
        z[0] = -x;
        z[t] = E::Fr::one();
        P::from_coefficients_vec(z)
    }

    // r(X) = f1(x) + f2(x)X + ... + ft(x)X^{t-1}
    fn fflonk_remainder_poly(fs: &[P], x: P::Point) -> P {
        let vs = fs.iter().map(|fi| fi.evaluate(&x)).collect();
        P::from_coefficients_vec(vs)
    }

    fn multi_eval(poly: &P, points: &[P::Point]) -> Vec<P::Point> {
        assert_eq!(poly.degree() + 1, points.len());
        points.iter().map(|p| poly.evaluate(p)).collect()
    }

    fn fflonk_aggregate_remainder(rxs: &[P::Point], vss: Vec<Vec<P::Point>>, t: usize) -> P {
        assert!(vss.iter().all(|vs| vs.len() == t));
        let rootss = rxs.iter().map(|&rx| Self::fflonk_roots(rx, t));
        let polys = vss.iter().map(|vs| P::from_coefficients_slice(vs));
        assert_eq!(rootss.len(), polys.len());
        let xs: Vec<_> = rootss.clone().flatten().collect();
        let ys: Vec<_> = polys.zip(rootss).flat_map(|(poly, roots)| Self::multi_eval(&poly, &roots)).collect();
        Self::interpolate(&xs, &ys)
    }

    fn verify_fflonk_1_poly_n_points(
        evk: &ExtendedVerifierKey<E>,
        g1: E::G1Affine,
        q1: E::G1Affine,
        t_root_xs: &[P::Point],
        vss: Vec<Vec<P::Point>>,
        t: usize,
    ) -> bool {
        let r = Self::fflonk_aggregate_remainder(t_root_xs, vss, t);
        let r1 = Self::commit_g(&evk.powers_in_g1, &r);

        let xs: Vec<_> = t_root_xs.iter().map(|&t_root_x| t_root_x.pow([t as u64])).collect();
        let z = Self::flonk_n_point_z(&xs, t);
        let z2 = Self::commit_g(&evk.powers_in_g2, &z);

        Self::verify(g1, r1, q1, z2, evk.powers_in_g2[0])
    }


    // Verifies single fflonk proof using n-point KZG.
    // Checks g - r = q * z in the exponent,
    // g = combine(f1,...,ft) is committed to G1,
    // r(X) = f1(x) + f2(x)X + ... + ft(x)X^{t-1}, // fflonk_remainder_poly
    // z(X) = X^t - x, // fflonk_vanishing_poly
    // and q = g // z = (g - r) / z is commited to G1. // fflonk_vanishing_poly
    fn verify_fflonk(
        evk: &ExtendedVerifierKey<E>,
        g1: E::G1Affine,
        q1: E::G1Affine,
        x: P::Point,
        vs: &[P::Point],
        t: usize,
    ) -> bool {
        let r = P::from_coefficients_slice(vs);
        let r1 = Self::commit_g(&evk.powers_in_g1, &r);
        let z = Self::fflonk_vanishing_poly(t, x);
        let z2 = Self::commit_g(&evk.powers_in_g2, &z);

        Self::verify(g1, r1, q1, z2, evk.powers_in_g2[0])
    }

    // checks f - r = q * z in the exponent given q and f are committed to G1
    fn verify(
        f1: E::G1Affine,
        r1: E::G1Affine,
        q1: E::G1Affine,
        z2: E::G2Affine,
        g2: E::G2Affine, // generator
    ) -> bool {
        let f1 = f1.into_projective();
        let r1 = r1.into_projective();
        let fr = f1 - r1;
        assert_eq!(f1, fr + r1);
        E::product_of_pairings(&[
            (fr.into_affine().into(), (-g2).into()),
            (q1.into(), z2.into()),
        ]).is_one()
    }

    fn interpolate(xs: &[E::Fr], ys: &[E::Fr]) -> P {
        let x1 = xs[0];
        let mut l = Self::z_of_point(x1);
        for &xj in xs.iter().skip(1) {
            let q = Self::z_of_point(xj);
            l = &l * &q;
        }

        let mut ws = vec![];
        for xj in xs {
            let mut wj = E::Fr::one();
            for xk in xs {
                if xk != xj {
                    let d = *xj - xk;
                    wj *= d;
                }
            }
            ws.push(wj);
        }
        ark_ff::batch_inversion(&mut ws);

        let mut res = P::zero();
        for ((&wi, &xi), &yi) in ws.iter().zip(xs).zip(ys) {
            let d = Self::z_of_point(xi);
            let mut z = &l / &d;
            z = &z * wi;
            z = &z * yi;
            res = res + z;
        }
        res
    }

    fn divide_with_remainder(a: &P, b: &P) -> (P, P) {
        let a: DenseOrSparsePolynomial<E::Fr> = a.into();
        let b: DenseOrSparsePolynomial<E::Fr> = b.into();
        let (q, r) = a.divide_with_q_and_r(&b).unwrap();
        (q.into(), r.into())
    }
}

fn skip_leading_zeros_and_convert_to_bigints<F: PrimeField, P: UVPolynomial<F>>(
    p: &P,
) -> (usize, Vec<F::BigInt>) {
    let mut num_leading_zeros = 0;
    while num_leading_zeros < p.coeffs().len() && p.coeffs()[num_leading_zeros].is_zero() {
        num_leading_zeros += 1;
    }
    let coeffs = convert_to_bigints(&p.coeffs()[num_leading_zeros..]);
    (num_leading_zeros, coeffs)
}

fn convert_to_bigints<F: PrimeField>(p: &[F]) -> Vec<F::BigInt> {
    let to_bigint_time = start_timer!(|| "Converting polynomial coeffs to bigints");
    let coeffs = ark_std::cfg_iter!(p)
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    end_timer!(to_bigint_time);
    coeffs
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;
    use ark_poly::Polynomial;
    use ark_poly::univariate::DensePolynomial;
    use ark_bw6_761::{BW6_761, Fr};

    use rand::Rng;
    use ark_ff::Field;
    use ark_ec::group::Group;
    use std::iter::FromIterator;


    type Bw6Poly = DensePolynomial<<BW6_761 as PairingEngine>::Fr>;
    type KzgBw6 = KZG10<BW6_761, Bw6Poly>;


    fn _setup_commit_open_check_test<E, P>(max_degree: usize)
        where
            E: PairingEngine,
            P: UVPolynomial<E::Fr, Point=E::Fr> + From<DensePolynomial<E::Fr>>,
            for<'a, 'b> &'a P: Div<&'b P, Output=P>,
            for<'a, 'b> &'a P: Mul<&'b P, Output=P>,
            for<'a, 'b> &'a P: Sub<&'b P, Output=P>,
            for<'a> &'a P: Mul<E::Fr, Output=P>,
            for<'a> &'a P: Into<DenseOrSparsePolynomial<'a, E::Fr>>,
    {
        let rng = &mut test_rng();

        let params = KZG10::<E, P>::setup(max_degree, rng);
        let pk = params.get_pk();
        let pvk = params.get_vk().prepare();

        let degree = rng.gen_range(0..=max_degree);
        let poly = P::rand(degree, rng);
        let comm = KZG10::<E, P>::commit(&pk, &poly);

        let x = E::Fr::rand(rng);
        let proof = KZG10::<E, P>::open(&pk, &poly, x);
        let y = poly.evaluate(&x);
        let valid = KZG10::<E, P>::check(&pvk, &comm, x, y, proof);
        assert!(valid, "check failed for a valid point");
        let not_y = E::Fr::rand(rng);
        let valid = KZG10::<E, P>::check(&pvk, &comm, x, not_y, proof);
        assert!(!valid, "check failed for an invalid point");
    }

    #[test]
    fn setup_commit_open_check_test() {
        _setup_commit_open_check_test::<BW6_761, Bw6Poly>(0);
        _setup_commit_open_check_test::<BW6_761, Bw6Poly>(1);
        _setup_commit_open_check_test::<BW6_761, Bw6Poly>(123);
    }

    #[test]
    fn test_fflonk_vanishing_poly() {
        let rng = &mut test_rng();

        let t = 2;
        let z = Fr::rand(rng);
        let x = z.pow([t as u64]);
        let vp = KzgBw6::fflonk_vanishing_poly(t, x);
        let zws = KzgBw6::fflonk_roots(z, t);
        let vp2 = zws.iter()
            .map(|&zw| KzgBw6::z_of_point(zw))
            .reduce(|p, pi| &p * &pi)
            .unwrap();
        assert_eq!(vp, vp2);
    }

    #[test]
    fn test_fflonk() {
        let rng = &mut test_rng();

        let t: usize = 4;
        let d = 15;
        let z = Fr::rand(rng);
        let x = z.pow([t as u64]);

        let fs: Vec<_> = (0..t)
            .map(|_| DensePolynomial::rand(d, rng))
            .collect();

        let g = KzgBw6::combine(&fs);
        let q = KzgBw6::fflonk_quotient_poly(&g, t, x);
        let r = KzgBw6::fflonk_remainder_poly(&fs, x);
        let z = KzgBw6::fflonk_vanishing_poly(t, x);

        let p = &g - &r;
        let (q2, r) = KzgBw6::divide_with_remainder(&p, &z);
        assert!(r.is_zero());
        assert_eq!(q, q2);
    }

    #[test]
    fn test_shplonk2_n_points() {
        let rng = &mut test_rng();
        let max_degree = 10;
        let k = 4;

        let params = KzgBw6::setup(max_degree, rng);
        let pk = params.get_pk();
        let pvk = params.get_vk().prepare();

        let degree = rng.gen_range(0..=max_degree);
        let poly = DensePolynomial::rand(degree, rng);
        let comm = KzgBw6::commit(&pk, &poly);

        let xs: Vec<_> = (0..k).map(|_| Fr::rand(rng)).collect();

        let xs_set = HashSet::from_iter(xs.clone());
        let (qf, ql) = KzgBw6::shplonk_open_1_poly_n_points_interactive(&pk, &poly, &xs_set, &mut Transcript::new(b"shplonk-test"));

        let ys: Vec<_> = xs.iter().map(|x| poly.evaluate(x)).collect();
        let valid = KzgBw6::shplonk_verify_1_poly_n_points_interactive(&pvk, &comm, qf, ql, &xs, &ys, &mut Transcript::new(b"shplonk-test"));
        assert!(valid, "check failed for a valid point");
    }

    #[test]
    fn test_shplonk() {
        let rng = &mut test_rng();

        let max_degree = 10;
        let n = 4;

        let params = KzgBw6::setup(max_degree, rng);
        let pk = params.get_pk();
        let pvk = params.get_vk().prepare();

        let fs: Vec<_> = (0..n).map(|_| DensePolynomial::rand(max_degree, rng)).collect();
        let f1s: Vec<_> = fs.iter().map(|fi| KzgBw6::commit(&pk, &fi)).collect();

        let (xss, yss): (Vec<_>, Vec<_>) = fs.iter().map(|fi| {
            let xs = (0..2).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
            let ys = xs.iter().map(|xi| fi.evaluate(&xi)).collect::<Vec<_>>();
            (xs, ys)
        }).unzip();

        let sets: Vec<_> = xss.iter().cloned().map(|xs| HashSet::from_iter(xs)).collect();

        let (qf, ql) = KzgBw6::shplonk_open(&pk, &fs, &sets, &mut Transcript::new(b"shplonk-test"));

        let valid = KzgBw6::shplonk_verify(&pvk, &f1s, qf, ql, &xss, &yss, &mut Transcript::new(b"shplonk-test"));
        assert!(valid);
    }

    #[test]
    fn test_fflonk_with_n_point_kzg() {
        let rng = &mut test_rng();

        let t: usize = 4;
        let d = 15;

        let max_degree = t * (d + 1);

        let params = KzgBw6::setup(max_degree, rng);
        let pk = params.get_pk();
        let evk = params.get_evk();

        let fs: Vec<_> = (0..t)
            .map(|_| DensePolynomial::rand(d, rng))
            .collect();
        let g = KzgBw6::combine(&fs);

        let g1 = KzgBw6::commit(&pk, &g);

        let z = Fr::rand(rng);
        let x = z.pow([t as u64]);

        let q1 = KzgBw6::open_combined(&pk, &g, x, t);

        let mut vs: Vec<_> = fs.iter().map(|fi| fi.evaluate(&x)).collect();

        let valid = KzgBw6::verify_fflonk(&evk, g1, q1, x, &vs, t);
        assert!(valid, "check failed for a valid point");

        vs[0] = Fr::rand(rng);
        let valid = KzgBw6::verify_fflonk(&evk, g1, q1, x, &vs, t);
        assert!(!valid, "check failed for an invalid point");
    }

    #[test]
    fn test_fflonk_in_m_points_with_n_point_kzg() {
        let rng = &mut test_rng();

        let t: usize = 4;
        let d = 15;
        let n = 2;

        let max_degree = t * (d + 1);

        let params = KzgBw6::setup(max_degree, rng);
        let pk = params.get_pk();
        let evk = params.get_evk();

        let fs: Vec<_> = (0..t)
            .map(|_| DensePolynomial::rand(d, rng))
            .collect();
        let g = KzgBw6::combine(&fs);

        let g1 = KzgBw6::commit(&pk, &g);

        let zs: Vec<_> = (0..n).map(|_| Fr::rand(rng)).collect();
        let xs: Vec<_> = zs.iter().map(|z| z.pow([t as u64])).collect();

        let q1 = KzgBw6::flonk_open_1_poly_n_points(&pk, &g, &xs, t);

        let vss = xs.iter().map(|x| fs.iter().map(|fi| fi.evaluate(&x)).collect()).collect();

        let valid = KzgBw6::verify_fflonk_1_poly_n_points(&evk, g1, q1, &zs, vss, t);
        assert!(valid, "check failed for a valid point");
    }

    #[test]
    fn commitment_linearity_test()
    {
        let rng = &mut test_rng();
        let max_degree = 123;

        let params = KzgBw6::setup(max_degree, rng);
        let pk = params.get_pk();

        let poly1 = Bw6Poly::rand(rng.gen_range(0..=max_degree), rng);
        let poly2 = Bw6Poly::rand(rng.gen_range(0..=max_degree), rng);

        let comm1 = KzgBw6::commit(&pk, &poly1);
        let comm2 = KzgBw6::commit(&pk, &poly2);

        let r = Fr::rand(rng);
        let agg_poly = KzgBw6::aggregate_polynomials(r, &[poly1, poly2]);
        let agg_comm = KzgBw6::aggregate_commitments(r, &[comm1, comm2]);

        let agg_poly_comm = KzgBw6::commit(&pk, &agg_poly);
        assert_eq!(agg_comm, agg_poly_comm, "Commit(af+bg) == aCommit(f) + bCommit(g)");
    }

    #[test]
    fn additivity_test()
    {
        let rng = &mut test_rng();
        let max_degree = 123;

        let params = KzgBw6::setup(max_degree, rng);
        let pk = params.get_pk();
        let pvk = params.get_vk().prepare();

        let poly1 = Bw6Poly::rand(rng.gen_range(0..=max_degree), rng);
        let poly2 = Bw6Poly::rand(rng.gen_range(0..=max_degree), rng);
        let x = Fr::rand(rng);
        let y1 = poly1.evaluate(&x);
        let y2 = poly2.evaluate(&x);

        // prover and verifier both know the randomizer
        let r = Fr::rand(rng);

        // prover
        let comm1 = KzgBw6::commit(&pk, &poly1);
        let comm2 = KzgBw6::commit(&pk, &poly2);
        let agg_poly = KzgBw6::aggregate_polynomials(r, &[poly1.clone(), poly2.clone()]);
        let agg_proof = KzgBw6::open(&pk, &agg_poly, x);

        // verifier
        let agg_value = KzgBw6::aggregate_values(r, &[y1, y2]);
        let agg_comm = KzgBw6::aggregate_commitments(r, &[comm1, comm2]);
        let valid = KzgBw6::check(&pvk, &agg_comm, x, agg_value, agg_proof);
        assert!(valid, "check failed for a valid point");

        // and also // TODO: a separate test
        let proof1 = KzgBw6::open(&pk, &poly1, x);
        let proof2 = KzgBw6::open(&pk, &poly2, x);
        assert_eq!(agg_proof, proof1 + proof2.mul(r).into_affine());
    }

//
//
//
//     for _ in 0..100 {
//         let mut degree = 0;
//         while degree <= 1 {
//             degree = usize::rand(rng) % 20;
//         }
//         let pp = KZG10::<E, P>::setup(degree, false, rng)?;
//         let (ck, vk) = KZG10::<E, P>::trim(&pp, degree)?;
//         let p = P::rand(degree, rng);
//         let hiding_bound = Some(1);
//         let (comm, rand) = KZG10::<E, P>::commit(&ck, &p, hiding_bound, Some(rng))?;
//         let point = E::Fr::rand(rng);
//         let value = p.evaluate(&point);
//         let proof = KZG10::<E, P>::open(&ck, &p, point, &rand)?;
//         assert!(
//             KZG10::<E, P>::check(&vk, &comm, point, value, &proof)?,
//             "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
//             degree,
//             p.degree(),
//             hiding_bound,
//         );
//     }
//     Ok(())
// }
//
// impl<E: PairingEngine, P: UVPolynomial<E::Fr>> KZG10<E, P> {
//     /// Specializes the public parameters for a given maximum degree `d` for polynomials
//     /// `d` should be less that `pp.max_degree()`.
//     pub(crate) fn trim(
//         pp: &UniversalParams<E>,
//         mut supported_degree: usize,
//     ) -> Result<(Powers<E>, VerifierKey<E>), Error> {
//         if supported_degree == 1 {
//             supported_degree += 1;
//         }
//         let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
//         let powers_of_gamma_g = (0..=supported_degree)
//             .map(|i| pp.powers_of_gamma_g[&i])
//             .collect();
//
//         let powers = Powers {
//             powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
//             powers_of_gamma_g: ark_std::borrow::Cow::Owned(powers_of_gamma_g),
//         };
//         let vk = VerifierKey {
//             g: pp.powers_of_g[0],
//             gamma_g: pp.powers_of_gamma_g[&0],
//             h: pp.h,
//             beta_h: pp.beta_h,
//             prepared_h: pp.prepared_h.clone(),
//             prepared_beta_h: pp.prepared_beta_h.clone(),
//         };
//         Ok((powers, vk))
//     }
// }
//
// #[test]
// fn add_commitments_test() {
//     let rng = &mut test_rng();
//     let p = DensePoly::from_coefficients_slice(&[
//         Fr::rand(rng),
//         Fr::rand(rng),
//         Fr::rand(rng),
//         Fr::rand(rng),
//         Fr::rand(rng),
//     ]);
//     let f = Fr::rand(rng);
//     let mut f_p = DensePoly::zero();
//     f_p += (f, &p);
//
//     let degree = 4;
//     let pp = KZG_Bls12_381::setup(degree, false, rng).unwrap();
//     let (powers, _) = KZG_Bls12_381::trim(&pp, degree).unwrap();
//
//     let hiding_bound = None;
//     let (comm, _) = KZG10::commit(&powers, &p, hiding_bound, Some(rng)).unwrap();
//     let (f_comm, _) = KZG10::commit(&powers, &f_p, hiding_bound, Some(rng)).unwrap();
//     let mut f_comm_2 = Commitment::empty();
//     f_comm_2 += (f, &comm);
//
//     assert_eq!(f_comm, f_comm_2);
// }
//
// fn end_to_end_test_template<E, P>() -> Result<(), Error>
//     where
//         E: PairingEngine,
//         P: UVPolynomial<E::Fr, Point = E::Fr>,
//         for<'a, 'b> &'a P: Div<&'b P, Output = P>,
// {
//     let rng = &mut test_rng();
//     for _ in 0..100 {
//         let mut degree = 0;
//         while degree <= 1 {
//             degree = usize::rand(rng) % 20;
//         }
//         let pp = KZG10::<E, P>::setup(degree, false, rng)?;
//         let (ck, vk) = KZG10::<E, P>::trim(&pp, degree)?;
//         let p = P::rand(degree, rng);
//         let hiding_bound = Some(1);
//         let (comm, rand) = KZG10::<E, P>::commit(&ck, &p, hiding_bound, Some(rng))?;
//         let point = E::Fr::rand(rng);
//         let value = p.evaluate(&point);
//         let proof = KZG10::<E, P>::open(&ck, &p, point, &rand)?;
//         assert!(
//             KZG10::<E, P>::check(&vk, &comm, point, value, &proof)?,
//             "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
//             degree,
//             p.degree(),
//             hiding_bound,
//         );
//     }
//     Ok(())
// }
//
// fn linear_polynomial_test_template<E, P>() -> Result<(), Error>
//     where
//         E: PairingEngine,
//         P: UVPolynomial<E::Fr, Point = E::Fr>,
//         for<'a, 'b> &'a P: Div<&'b P, Output = P>,
// {
//     let rng = &mut test_rng();
//     for _ in 0..100 {
//         let degree = 50;
//         let pp = KZG10::<E, P>::setup(degree, false, rng)?;
//         let (ck, vk) = KZG10::<E, P>::trim(&pp, 2)?;
//         let p = P::rand(1, rng);
//         let hiding_bound = Some(1);
//         let (comm, rand) = KZG10::<E, P>::commit(&ck, &p, hiding_bound, Some(rng))?;
//         let point = E::Fr::rand(rng);
//         let value = p.evaluate(&point);
//         let proof = KZG10::<E, P>::open(&ck, &p, point, &rand)?;
//         assert!(
//             KZG10::<E, P>::check(&vk, &comm, point, value, &proof)?,
//             "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
//             degree,
//             p.degree(),
//             hiding_bound,
//         );
//     }
//     Ok(())
// }
//
// fn batch_check_test_template<E, P>() -> Result<(), Error>
//     where
//         E: PairingEngine,
//         P: UVPolynomial<E::Fr, Point = E::Fr>,
//         for<'a, 'b> &'a P: Div<&'b P, Output = P>,
// {
//     let rng = &mut test_rng();
//     for _ in 0..10 {
//         let mut degree = 0;
//         while degree <= 1 {
//             degree = usize::rand(rng) % 20;
//         }
//         let pp = KZG10::<E, P>::setup(degree, false, rng)?;
//         let (ck, vk) = KZG10::<E, P>::trim(&pp, degree)?;
//         let mut comms = Vec::new();
//         let mut values = Vec::new();
//         let mut points = Vec::new();
//         let mut proofs = Vec::new();
//         for _ in 0..10 {
//             let p = P::rand(degree, rng);
//             let hiding_bound = Some(1);
//             let (comm, rand) = KZG10::<E, P>::commit(&ck, &p, hiding_bound, Some(rng))?;
//             let point = E::Fr::rand(rng);
//             let value = p.evaluate(&point);
//             let proof = KZG10::<E, P>::open(&ck, &p, point, &rand)?;
//
//             assert!(KZG10::<E, P>::check(&vk, &comm, point, value, &proof)?);
//             comms.push(comm);
//             values.push(value);
//             points.push(point);
//             proofs.push(proof);
//         }
//         assert!(KZG10::<E, P>::batch_check(
//             &vk, &comms, &points, &values, &proofs, rng
//         )?);
//     }
//     Ok(())
// }
//
// #[test]
// fn end_to_end_test() {
//     end_to_end_test_template::<Bls12_377, UniPoly_377>().expect("test failed for bls12-377");
//     end_to_end_test_template::<Bls12_381, UniPoly_381>().expect("test failed for bls12-381");
// }
//
// #[test]
// fn linear_polynomial_test() {
//     linear_polynomial_test_template::<Bls12_377, UniPoly_377>()
//         .expect("test failed for bls12-377");
//     linear_polynomial_test_template::<Bls12_381, UniPoly_381>()
//         .expect("test failed for bls12-381");
// }
// #[test]
// fn batch_check_test() {
//     batch_check_test_template::<Bls12_377, UniPoly_377>().expect("test failed for bls12-377");
//     batch_check_test_template::<Bls12_381, BW6_POLY>().expect("test failed for bls12-381");
// }
}
