use crate::piop::bitmask_packing::{BitmaskPackingRegisters, SuccinctAccountableRegisterEvaluations, BitmaskPackingPolynomials};
use crate::piop::Protocol;
use crate::domains::Domains;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_bls12_377::G1Affine;
use crate::Bitmask;
use ark_bw6_761::Fr;
use crate::piop::affine_addition::{AffineAdditionRegisters, PartialSumsPolynomials};

pub struct PackedRegisterBuilder {
    bitmask: Bitmask,
    affine_addition_registers: AffineAdditionRegisters,
    bitmask_packing_registers: Option<BitmaskPackingRegisters>,
    register_evaluations: Option<SuccinctAccountableRegisterEvaluations>,
}

impl Protocol for PackedRegisterBuilder {
    type P1 = PartialSumsPolynomials;
    type P2 = BitmaskPackingPolynomials;
    type E = SuccinctAccountableRegisterEvaluations;

    fn init(domains: Domains, bitmask: Bitmask, pks: Vec<G1Affine>) -> Self {
        PackedRegisterBuilder {
            bitmask: bitmask.clone(),
            affine_addition_registers: AffineAdditionRegisters::new(domains, &bitmask, pks),
            bitmask_packing_registers: None,
            register_evaluations: None,
        }
    }

    fn get_1st_round_register_polynomials(&self) -> PartialSumsPolynomials {
        self.affine_addition_registers.get_partial_sums_register_polynomials()
    }

    fn get_2nd_round_register_polynomials(&mut self, bitmask_chunks_aggregation_challenge: Fr) -> BitmaskPackingPolynomials {
        let bitmask_packing_registers = BitmaskPackingRegisters::new(
            self.affine_addition_registers.domains.clone(),
            &self.bitmask,
            bitmask_chunks_aggregation_challenge
        );
        let res = bitmask_packing_registers.get_register_polynomials();
        self.bitmask_packing_registers = Some(bitmask_packing_registers);
        res
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<Fr>> {
        let affine_addition_constraints = self.affine_addition_registers.compute_constraint_polynomials();
        let bitmask_packing_constraints = self.bitmask_packing_registers.as_ref().unwrap().compute_constraint_polynomials();
        let mut constraints = vec![];
        constraints.extend(affine_addition_constraints);
        constraints.extend(bitmask_packing_constraints);
        constraints
    }

    fn evaluate_register_polynomials(&mut self, point: Fr) -> SuccinctAccountableRegisterEvaluations {
        let affine_addition_evals = self.affine_addition_registers.evaluate_register_polynomials(point);
        let bitmask_packing_evals = self.bitmask_packing_registers.as_ref().unwrap().evaluate_register_polynomials(point);
        let evals = SuccinctAccountableRegisterEvaluations {
            c: bitmask_packing_evals.0,
            acc: bitmask_packing_evals.1,
            basic_evaluations: affine_addition_evals,
        };
        self.register_evaluations = Some(evals.clone());
        evals
    }

    fn compute_linearization_polynomial(&self, phi: Fr, zeta_minus_omega_inv: Fr) -> DensePolynomial<Fr> {
        let affine_addition_lp =
            self.affine_addition_registers.compute_linearization_polynomial(&self.register_evaluations.as_ref().unwrap().basic_evaluations, phi, zeta_minus_omega_inv);
        let bitmask_packing_lp =
            self.bitmask_packing_registers.as_ref().unwrap().compute_linearization_polynomial(self.register_evaluations.as_ref().unwrap(), phi, zeta_minus_omega_inv);
        affine_addition_lp + bitmask_packing_lp
    }

    fn get_all_register_polynomials(self) -> Vec<DensePolynomial<Fr>> {
        let affine_addition_polys = self.affine_addition_registers.get_register_polynomials().to_vec();
        let bitmask_packing_polys = self.bitmask_packing_registers.unwrap().get_register_polynomials().to_vec();
        let mut polys = vec![];
        polys.extend(affine_addition_polys);
        polys.extend(bitmask_packing_polys);
        polys
    }
}