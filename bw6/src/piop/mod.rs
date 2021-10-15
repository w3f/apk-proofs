use ark_bw6_761::{Fr, G1Affine};
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use crate::{Bitmask, utils, PublicInput, Keyset};
use ark_poly::Radix2EvaluationDomain;
use crate::domains::Domains;


pub mod affine_addition;
pub mod bitmask_packing;
pub mod bit_counting;

pub mod basic;
pub mod packed;
pub mod counting;


pub trait RegisterCommitments: CanonicalSerialize + CanonicalDeserialize {
    fn as_vec(&self) -> Vec<G1Affine>;
}

pub trait RegisterPolynomials {
    type C: RegisterCommitments;
    fn commit<F: Clone + Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, f: F) -> Self::C;
}

impl RegisterCommitments for () {
    fn as_vec(&self) -> Vec<G1Affine> {
        vec![]
    }
}

impl RegisterPolynomials for () {
    type C = ();

    fn commit<F: Fn(&DensePolynomial<Fr>) -> G1Affine>(&self, _f: F) -> Self::C {
        ()
    }
}

// PIOP stays for Polynomial Interactive Oracle Proof.
// It's a type of information-theoretical interactive proof systems where to convince the verifier that
// a relation holds, the prover sends bounded-degree polynomials, that are accessible by the verifier
// as "oracles", meaning it can query them for evaluations in some public-coin points.
// A PIOP can be lifted to a SNARK by implementing the oracle queries with a polynomial
// commitment scheme under some computational assumption and performing a Fiat-Shamir transform.

// ProverProtocol and VerifierProtocol together encapsulate a PIOP optimized for a limited class of computations.
// The computation is expressed in a polynomial form in a way similar to AIR arithmetization:

// 0. Following PLONK we fix a "domain", a subgroup in the multiplicative group of the field of size n and its generator w.
//    Now we can interpret an array of length n as a degree < n polynomial, by interpolating over the domain:
//    [p_0,...,p_{n-1}] is equivalent to the unique polynomial p of degree < n, such that p(w^i) = p_i, i = 0,...,n-1.
//    Polynomial p'(Z) = p(Zw) represents the left circular shift [p_1,...,p_0] of the array:
//    p'(w^i) = p(w^{i+1}) = p_{i+1}, i = 0,...,n-2, and p'(w^{n-1}) = p(w^n) = p(w^0) = p_0.

// 1. The computation state is represented as registers, arrays of length n. Together registers can be visualized as a table,
//    where each row corresponds to the computation state at some step. Registers give raise to register polynomials.
//    The consecutive element of a register r can be addressed using the shifted polynomial r(Zw).

// 2. The computation is defined by a constraint polynomial f. It is a univariate polynomial of the form:
//    f(Z) = P(r1(Z),...,rm(Z),r1(Zw),...,rm(Zw)), where P is a small-degree multivariate polynomial that is known to the verifier.
//    degree of f is then <= (n-1)deg(P)
//    The computation is valid iff the constraint polynomial is zero for every 2 consecutive rows:
//    f(w^i) = f(r1(w^i),...,rm(w^i),r1(w^{i+1}),...,rm(w^i{i+1})) = 0, i = 0,...,n-1.
//    Notice that the constraint should also hold between the last and the first rows of the state table.

// The formula above says that the computation is valid iff f = 0 over the domain or, in other words, f has a root at every point of the domain.
// Recall that z_0 is a root of f iff the polynomial (Z-z_0) divides f(Z), i.e. f(Z) = q(Z)(Z-z_0). Inductively we get that
// the computation is valid iff there exists a polynomial q such that f(Z) = q(Z)(Z-w^0)...(Z-w^{n-1}) = q(Z)(Z^n - 1).

// A PIOP for the arithmetization described above is:
// 0. Verifier has an oracle access to some "public input" registers, and knows P, a formula for the constraint polynomial.
//    Prover knows the same, and has a satisfying assignment to other registers,
//    consistent both with the public input registers and the constraint polynomial.
// 1. Prover sends register polynomials r1,...,rl, other than "public input" registers.
// 2. Prover computes the constraint polynomial f, the quotient polynomial q(Z) = f(Z) / (Z^n - 1) and sends q.
// 3. Verifier for a random point z:
//    - computes f(z) using the formula P and queries to the register polynomials in z and zw
//    - computes z^n - 1
//    - queries q in z and checks f(z) = q(z) * (z^n - 1) holds

// As the constraint polynomial degree is bounded, the soundness of the protocol follows from Schwartz–Zippel lemma.

// The implementation differs from the abstract protocol described above in the following.
// 1. It's essential to have multiple constraints f_i such that each f_i is 0. To reduce this case to the protocol above,
//    there's an additional round of interaction: after the prover sends the registers, it receives a random challenge phi
//    from the verifier, that is used to aggregate the constraints: f = \sum phi^i f_i
// 2. There's another additional round of interaction specific to the "packed" case.
// 3. To reduce communication (oracle queries), instead of querying multiple registers in zw,
//    verifier is given the "linearization" polynomial that is enough to be queried once.
//    It is efficient if the constraint polynomial is linear in all "shifted" terms ri(Zw).

pub trait ProverProtocol {
    type P1: RegisterPolynomials;
    type P2: RegisterPolynomials;
    type E: RegisterEvaluations;
    type PI: PublicInput;

    fn init(domains: Domains, bitmask: Bitmask, keyset: Keyset) -> Self;

    // These 2 methods together return register polynomials the prover should commit to.
    // The 2nd one is used only in the "packed" scheme as it requires an additional challenge
    // (to aggregate the bitmask chunks) from the verifier,
    // that can be received only after the bitmask has been committed.
    fn get_register_polynomials_to_commit1(&self) -> Self::P1;
    fn get_register_polynomials_to_commit2(&mut self, verifier_challenge: Fr) -> Self::P2;

    // This method returns register polynomials the prover should open. Those are the same polynomials
    // as the previous 2 methods together, and additionally 2 polynomials representing the keyset
    // (prover doesn't need to commit to them, as verifier knows them anyway, but still should open).
    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<Fr>>;

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>>;

    //TODO: remove domains param
    fn compute_quotient_polynomial(&self, phi: Fr, domain: Radix2EvaluationDomain<Fr>) -> DensePolynomial<Fr> {
        let w = utils::randomize(phi, &self.compute_constraint_polynomials());
        let (q_poly, r) = w.divide_by_vanishing_poly(domain).unwrap();
        assert_eq!(r, DensePolynomial::zero());
        q_poly
    }

    fn evaluate_register_polynomials(&mut self, point: Fr) -> Self::E;

    // Some constraints require access to 2 consecutive rows of the registers, e.g. affine addition
    // adds 2 points that are represented as consecutive rows. In polynomials it is represented by evaluating
    // a register polynomial in 2 points: zeta and zeta * omega. (An alternative would be for each register f(Z),
    // that needs to be constrained in more than one location, introduce an additional "shifted" register f(Z * omega),
    // and evaluate all the registers in the single point. But proving evaluations is such way would require additionally
    // committing to the new registers.) So the prover batch-opens all the registers in zeta, but instead of also batch-openning some of them in zeta * omega
    // (that would require communicating more evaluations), it opens in zeta * omega only their linear combination.
    // Example: //TODO
    // The verifier can restore the commitment to this "linearization" polynomial from the commitments to the register polynomials and their evaluations in zeta,
    // so the required communication (for any number of polynomials) is just the proof and the evaluation.
    // Plonk section "Reducing the number of field elements" describes the same for some more general case.
    fn compute_linearization_polynomial(&self, phi: Fr, zeta: Fr) -> DensePolynomial<Fr>;
}

pub trait RegisterEvaluations: CanonicalSerialize + CanonicalDeserialize {
    fn as_vec(&self) -> Vec<Fr>;
}

pub trait VerifierProtocol {
    type C1: RegisterCommitments; // commitments to ProverProtocol::P1
    type C2: RegisterCommitments; // commitments to ProverProtocol::P2

    fn restore_commitment_to_linearization_polynomial(
        &self,
        phi: Fr,
        zeta_minus_omega_inv: Fr,
        commitments: &Self::C1,
        extra_commitments: &Self::C2,
    ) -> ark_bw6_761::G1Projective;

    // fn evaluate_constraint_polynomials(
    //     &self,
    //     apk: ark_bls12_377::G1Affine,
    //     evals_at_zeta: &LagrangeEvaluations<Fr>,
    //     r: Fr,
    //     bitmask: &Bitmask,
    //     domain_size: u64,
    // ) -> Vec<Fr>;
}
