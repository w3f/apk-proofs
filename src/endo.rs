use ark_ff::{field_new, Zero, BitIteratorBE};
use ark_ff::biginteger::{BigInteger384, BigInteger768};
use ark_ec::ProjectiveCurve;
use ark_ec::bls12::Bls12Parameters;
use ark_bw6_761::{Fq, Fr, G1Affine, G1Projective};
use std::ops::AddAssign;

// See https://github.com/celo-org/zexe/blob/master/algebra/src/bw6_761/curves/g1.rs#L37-L71
// and also https://github.com/celo-org/zexe/blob/master/scripts/glv_lattice_basis/src/lib.rs

/// phi((x, y)) = (\omega x, y)
/// \omega = 0x531dc16c6ecd27aa846c61024e4cca6c1f31e53bd9603c2d17be416c5e44
/// 26ee4a737f73b6f952ab5e57926fa701848e0a235a0a398300c65759fc4518315
/// 1f2f082d4dcb5e37cb6290012d96f8819c547ba8a4000002f962140000000002a
const OMEGA: Fq = field_new!(
        Fq,
        BigInteger768([
            7467050525960156664,
            11327349735975181567,
            4886471689715601876,
            825788856423438757,
            532349992164519008,
            5190235139112556877,
            10134108925459365126,
            2188880696701890397,
            14832254987849135908,
            2933451070611009188,
            11385631952165834796,
            64130670718986244
        ])
    );

/// lambda in Z s.t. phi(P) = lambda*P for all P
/// \lambda = 0x9b3af05dd14f6ec619aaf7d34594aabc5ed1347970dec00452217cc900000008508c00000000001
const LAMBDA: Fr = field_new!(
        Fr,
        (BigInteger384([
            15766275933608376691,
            15635974902606112666,
            1934946774703877852,
            18129354943882397960,
            15437979634065614942,
            101285514078273488
        ]))
    );

const U: &'static [u64] = ark_bls12_377::Parameters::X;

fn mul_by_u(p: &G1Projective) -> G1Projective {
    let mut res = G1Projective::zero();
    for i in BitIteratorBE::without_leading_zeros(U) {
        res.double_in_place();
        if i {
            res.add_assign(p)
        }
    }
    res
}

fn glv_endomorphism_in_place(p: &mut G1Affine) {
    let x =  &mut p.x;
    *x *= &OMEGA;
}

fn glv_endomorphism(p: &G1Affine) -> G1Affine {
    G1Affine::new(p.x * OMEGA, p.y, false)
}

fn glv_endomorphism_proj(p: &G1Projective) -> G1Projective {
    G1Projective::new(p.x * OMEGA, p.y, p.z)
}

// See https://eprint.iacr.org/2020/351.pdf, section 3.1
pub fn subgroup_check(p: &G1Projective) -> bool {
    let up = mul_by_u(p);
    let u2p = mul_by_u(&up);
    let u3p = mul_by_u(&u2p);

    (up + p + glv_endomorphism_proj(&(u3p - u2p + p))).is_zero()
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{test_rng, UniformRand, Field, One};
    use ark_ec::AffineCurve;
    use ark_ec::models::short_weierstrass_jacobian::GroupAffine;

    #[test]
    pub fn test_omega() {
        assert!(OMEGA.pow([3]).is_one());
    }

    #[test]
    pub fn test_endo() {
        let rng = &mut test_rng();

        let p1 = ark_bw6_761::G1Projective::rand(rng).into_affine();
        let mut p2 = p1.clone();

        assert_eq!(glv_endomorphism(&p1), p1.mul(LAMBDA));
        glv_endomorphism_in_place(&mut p2);
        assert_eq!(p2, p1.mul(LAMBDA));
    }

    #[test]
    pub fn test_endo_proj() {
        let rng = &mut test_rng();

        let p = ark_bw6_761::G1Projective::rand(rng);

        assert_eq!(glv_endomorphism_proj(&p), p.mul(LAMBDA));
    }

    #[test]
    pub fn test_subgroup_check() {
        let rng = &mut test_rng();

        let p = ark_bw6_761::G1Projective::rand(rng);

        assert!(subgroup_check(&p));

        let point_not_in_g1 = loop {
            let x = Fq::rand(rng);
            let p = G1Affine::get_point_from_x(x, false);
            if p.is_some() && !p.unwrap().is_in_correct_subgroup_assuming_on_curve() {
                break p.unwrap();
            }
        };

        assert!(point_not_in_g1.is_on_curve());
        assert!(!point_not_in_g1.is_in_correct_subgroup_assuming_on_curve());
        assert!(!subgroup_check(&point_not_in_g1.into_projective()));
    }
}