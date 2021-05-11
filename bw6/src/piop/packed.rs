use crate::piop::bitmask_packing::{BitmaskPackingRegisters, SuccinctAccountableRegisterEvaluations, BitmaskPackingPolynomials};
use crate::piop::ProverProtocol;
use crate::domains::Domains;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_bls12_377::G1Affine;
use crate::{Bitmask, utils};
use ark_bw6_761::Fr;
use crate::piop::affine_addition::{AffineAdditionRegisters, PartialSumsPolynomials, PartialSumsAndBitmaskPolynomials};

pub struct PackedRegisterBuilder {
    bitmask: Bitmask,
    affine_addition_registers: AffineAdditionRegisters,
    bitmask_packing_registers: Option<BitmaskPackingRegisters>,
    register_evaluations: Option<SuccinctAccountableRegisterEvaluations>,
}

impl ProverProtocol for PackedRegisterBuilder {
    type P1 = PartialSumsAndBitmaskPolynomials;
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

    fn get_register_polynomials_to_commit1(&self) -> PartialSumsAndBitmaskPolynomials {
        self.affine_addition_registers.get_partial_sums_and_bitmask_polynomials()
    }


    fn get_register_polynomials_to_commit2(&mut self, bitmask_chunks_aggregation_challenge: Fr) -> BitmaskPackingPolynomials {
        let bitmask_packing_registers = BitmaskPackingRegisters::new(
            self.affine_addition_registers.domains.clone(),
            &self.bitmask,
            bitmask_chunks_aggregation_challenge,
        );
        let res = bitmask_packing_registers.get_register_polynomials();
        self.bitmask_packing_registers = Some(bitmask_packing_registers);
        res
    }

    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<Fr>> {
        let affine_addition_polys = self.affine_addition_registers.get_register_polynomials().to_vec();
        let bitmask_packing_polys = self.bitmask_packing_registers.unwrap().get_register_polynomials().to_vec();
        let mut polys = vec![];
        polys.extend(affine_addition_polys);
        polys.extend(bitmask_packing_polys);
        polys
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
        let evals = self.register_evaluations.as_ref().unwrap();
        let affine_addition_parts =
            self.affine_addition_registers.compute_constraints_linearized(&evals.basic_evaluations, zeta_minus_omega_inv);
        let bitmask_packing_parts =
            self.bitmask_packing_registers.as_ref().unwrap().compute_constraints_linearized(evals, zeta_minus_omega_inv);
        let mut parts = vec![];
        parts.extend(affine_addition_parts);
        parts.extend(bitmask_packing_parts);
        utils::randomize(phi, &parts)
    }
}