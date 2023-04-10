use curve25519_dalek::constants::{BASEPOINT_ORDER, RISTRETTO_BASEPOINT_TABLE};
use curve25519_dalek::ristretto::{CompressedRistretto, RistrettoPoint};
use curve25519_dalek::scalar::Scalar;
use once_cell::sync::Lazy;
use rand::{CryptoRng, RngCore};
use rug::{Complete, Integer, integer::Order, ops::RemRounding};
use sha3::{Shake128, digest::{ExtendableOutput, Update}};
use subtle::ConstantTimeEq;

use crate::gen_prime::*;

static CURVE_ORDER: Lazy<Integer> = Lazy::new(|| {
    Integer::from_digits(BASEPOINT_ORDER.as_bytes(), Order::Lsf)
});

// ceil(ln(2) * 2**129)
static ETA_SCALE_FACTOR: Lazy<Integer> = Lazy::new(|| {
    Integer::parse("471731526451026588275888285528308968799").unwrap().complete()
});

#[allow(non_snake_case)]
pub struct Prover {
    prf_key: [u8; 16],
    kP_pm1: Integer,
    kP_p: Integer,
    kP_qm1: Integer,
    kP_q: Integer,
    NP: Integer,
    NP_sq: Integer,
    NP_p: Integer,
    NP_p_sq: Integer,
    NP_q: Integer,
    NP_q_sq: Integer,
    NP_q_inv: Integer,
    NP_q_sq_inv: Integer,
    NV: Integer,
    g: Integer,
    g_inv: Integer,
    M: Integer,
}

#[allow(non_snake_case)]
pub struct Verifier {
    prf_key: [u8; 16],
    kV: Integer,
    NP: Integer,
    NP_sq: Integer,
    NV: Integer,
    g: Integer,
    Delta: Integer,
    Delta_scalar: Scalar,
    M_delta: Integer,
    u_hi_lim: Integer,
}

#[allow(non_snake_case)]
pub struct Proof {
    u_hi: Integer,
    s_inv: Integer, // Send (s^(-1)) instead of s, as it makes it easier for the verifier.
    R: CompressedRistretto,
    hash_t_V: [u8; 32],
}

#[derive(Debug)]
pub enum InvalidProof {
    InvalidCompressedPoint,
    OutOfRange,
    FailedConsistencyCheck,
}

impl Prover {
    #[allow(non_snake_case)]
    pub fn prove(&self, msg: &[u8]) -> (Scalar, RistrettoPoint, Proof) {
        let mut hasher = Shake128::default();
        hasher.update(b"0");
        hasher.update(msg);
        let base = hash_to_modulus(hasher.clone(), &self.NP_sq);

        hasher.update(&self.prf_key);
        let offset = hash_to_modulus(hasher, &self.NP);

        // Decompose into mod p and mod q, then use crt.
        let (base_p_add, base_p_pow_mul, base_p_pow_add) =
            pow_mod_p_sq(&base, &self.kP_pm1, &self.kP_p, &self.NP_p, &self.NP_p_sq);
        let (base_q_add, base_q_pow_mul, base_q_pow_add) =
            pow_mod_p_sq(&base, &self.kP_qm1, &self.kP_q, &self.NP_q, &self.NP_q_sq);
        let base_pow_mul = crt(base_p_pow_mul, base_q_pow_mul,
                               &self.NP_p_sq, &self.NP_q_sq, &self.NP_q_sq_inv);
        let base_pow_add = crt(base_p_pow_add, base_q_pow_add,
                               &self.NP_p, &self.NP_q, &self.NP_q_inv);

        let u = mod_round(&crt(base_p_add, base_q_add, &self.NP_p, &self.NP_q, &self.NP_q_inv), &self.NP);
        let v = mod_round(&(dist(base_pow_mul, &self.NP) + base_pow_add + offset), &self.NP);
        let (u_hi, u_lo) = u.div_rem_round_ref(&self.M).complete();

        let s_inv = self.g_inv.secure_pow_mod_ref(&u_lo, &self.NV).complete();
        let t = mod_round(&self.g.secure_pow_mod_ref(&v, &self.NV).complete(), &self.NV).abs();

        let r = integer_to_scalar(u_lo);
        let R = &r * RISTRETTO_BASEPOINT_TABLE;
        let V = &integer_to_scalar(v) * RISTRETTO_BASEPOINT_TABLE;

        hasher = Shake128::default();
        hasher.update(b"1");
        hash_integer(&mut hasher, &t, &self.NV);
        hasher.update(V.compress().as_bytes());

        let mut hash_t_V = [0u8; 32];
        hasher.finalize_xof_into(&mut hash_t_V);

        (r, R, Proof { u_hi, s_inv, R: R.compress(), hash_t_V })
    }
}

impl Verifier {
    #[allow(non_snake_case)]
    pub fn verify(&self, msg: &[u8], proof: Proof) -> Result<RistrettoPoint, InvalidProof> {
        let mut hasher = Shake128::default();
        hasher.update(b"0");
        hasher.update(msg);
        let base = hash_to_modulus(hasher.clone(), &self.NP_sq);

        hasher.update(&self.prf_key);
        let offset = hash_to_modulus(hasher, &self.NP);

        let w = dist(base.secure_pow_mod(&self.kV, &self.NP_sq), &self.NP) + offset;
        let w_lo = mod_round(&(w - &self.M_delta * &proof.u_hi), &self.NP);
        if proof.u_hi.abs() > self.u_hi_lim {
            return Err(InvalidProof::OutOfRange);
        }

        let R = proof.R.decompress().ok_or(InvalidProof::InvalidCompressedPoint)?;

        let t_check = mod_round(
            &(self.g.secure_pow_mod_ref(&w_lo, &self.NV).complete() *
            proof.s_inv.secure_pow_mod(&self.Delta, &self.NV)), &self.NV).abs();
        let V_check = &integer_to_scalar(w_lo) * RISTRETTO_BASEPOINT_TABLE - self.Delta_scalar * R;

        hasher = Shake128::default();
        hasher.update(b"1");
        hash_integer(&mut hasher, &t_check, &self.NV);
        hasher.update(V_check.compress().as_bytes());

        let mut hash_check = [0u8; 32];
        hasher.finalize_xof_into(&mut hash_check);
        if hash_check.ct_ne(&proof.hash_t_V).into() {
            return Err(InvalidProof::FailedConsistencyCheck);
        }

        Ok(R)
    }
}

#[allow(non_snake_case)]
pub fn setup<R: RngCore + CryptoRng>(rng: &mut R, ell: usize) -> (Prover, Verifier) {
    // Choose parameters.
    let (ell_prime, eta) = {
        let mut last_eta = Integer::from(1);
        loop {
            let ell_prime = (&last_eta * &*CURVE_ORDER).complete().significant_bits() as usize;
            let ell_prime = ell_prime + ell + 128 + 2;
            let ell_prime = ell_prime + (ell_prime % 2);

            let eta = &*ETA_SCALE_FACTOR * Integer::from(ell_prime);

            if eta <= last_eta {
                break (ell_prime, eta);
            } else {
                last_eta = eta;
            }
        }
    };

    assert!(eta < (Integer::from(1) << (ell / 2)));

    let mut prf_key = [0u8; 16];
    rng.fill_bytes(&mut prf_key);

    let NP_p = gen_prime(rng, ell_prime / 2);
    let NP_q = gen_prime(rng, ell_prime / 2);

    let NP_p_sq = (&NP_p * &NP_p).complete();
    let NP_q_sq = (&NP_q * &NP_q).complete();
    let NP = (&NP_p * &NP_q).complete();
    let NP_sq = (&NP_p_sq * &NP_q_sq).complete();

    let NP_q_inv = Integer::from(NP_q.invert_ref(&NP_p).unwrap());
    let NP_q_sq_inv = Integer::from(NP_q_sq.invert_ref(&NP_p_sq).unwrap());

    let phi: Integer = (&NP_p - Integer::from(1)) * (&NP_q - Integer::from(1));
    let d = &phi * Integer::from(phi.invert_ref(&NP).unwrap());

    let half_eta = &eta / Integer::from(2);
    let Delta = sample_mod(rng, &eta) - &half_eta;
    let kP = sample_mod(rng, &((&NP_sq * half_eta) << 128));
    let kV = &kP + d * &Delta;

    let NV_p = gen_safe_prime(rng, ell / 2);
    let NV_q = gen_safe_prime(rng, ell / 2);
    let NV = NV_p * NV_q;
    let g = sample_mod(rng, &NV);
    let g_inv = Integer::from(g.invert_ref(&NV).unwrap());
    let M = (&*CURVE_ORDER * &NV).complete();

    (
        Prover {
            prf_key,
            kP_pm1: &kP % (&NP_p - Integer::from(1)),
            kP_p: (&kP % &NP_p).complete(),
            kP_qm1: &kP % (&NP_q - Integer::from(1)),
            kP_q: (&kP % &NP_q).complete(),
            NP: NP.clone(),
            NP_sq: NP_sq.clone(),
            NP_p,
            NP_p_sq,
            NP_q,
            NP_q_sq,
            NP_q_inv,
            NP_q_sq_inv,
            NV: NV.clone(),
            g: g.clone(),
            g_inv,
            M: M.clone(),
        },
        Verifier {
            u_hi_lim: &NP / (Integer::from(2) * &M),
            M_delta: M * &Delta,
            prf_key,
            kV,
            NP,
            NP_sq,
            NV,
            g,
            Delta: Delta.clone(),
            Delta_scalar: integer_to_scalar(Delta),
        }
    )
}

// Decomposes x into 2 subgroups (multiplication mod p, addition mod p), then does the
// exponentiation inside each subgroup.
fn pow_mod_p_sq(x: &Integer, y_pm1: &Integer, y_p: &Integer, p: &Integer, p_sq: &Integer)
    -> (Integer, Integer, Integer) {
    let x = (x % p_sq).complete();
    let x_add_neg_exp = x.secure_pow_mod_ref(&(p - 1i8).complete(), p_sq).complete();
    let x_mul = (x * &x_add_neg_exp) % p_sq;
    let x_pow_mul = x_mul.secure_pow_mod(y_pm1, p_sq);

    let x_add = (1i8 - x_add_neg_exp).div_exact(p);
    let x_pow_add = (&x_add * y_p).complete();

    (x_add, x_pow_mul, x_pow_add)
}

fn crt(mod_p: Integer, mod_q: Integer, p: &Integer, q: &Integer, q_inv: &Integer) -> Integer {
    let h = (q_inv * (mod_p - &mod_q)).rem_floor(p);
    mod_q + q * h
}

fn hash_to_modulus<Hasher: ExtendableOutput>(hash: Hasher, modulus: &Integer) -> Integer {
    let num_bytes = ((modulus.significant_bits() + 128 + 7) / 8) as usize;
    Integer::from_digits(&hash.finalize_boxed(num_bytes), Order::Lsf) % modulus
}

fn sample_mod<R: RngCore + CryptoRng>(rng: &mut R, modulus: &Integer) -> Integer {
    let num_bytes = ((modulus.significant_bits() + 128 + 7) / 8) as usize;
    let mut bytes = vec![0u8; num_bytes];
    rng.fill_bytes(&mut bytes);
    Integer::from_digits(&bytes, Order::Lsf) % modulus
}

fn hash_integer<Hasher: Update>(hasher: &mut Hasher, x: &Integer, bound: &Integer) {
    let num_bytes = (bound.significant_bits() as usize + 7) / 8;
    let mut x_bytes = vec![0u8; num_bytes];
    x.write_digits(&mut x_bytes, Order::Lsf);
    hasher.update(&x_bytes);
}

fn dist(x: Integer, modulus: &Integer) -> Integer {
    let x_shifted = Integer::from(x.invert_ref(modulus).unwrap()) * x;
    (x_shifted - 1i8).div_exact(modulus)
}

fn mod_round(x: &Integer, modulus: &Integer) -> Integer {
    x.div_rem_round_ref(modulus).complete().1
}

fn integer_to_scalar(x: Integer) -> Scalar {
    let mut x_bytes = [0u8; 32];
    x.rem_floor(&*CURVE_ORDER).write_digits(&mut x_bytes, Order::Lsf);
    Scalar::from_canonical_bytes(x_bytes).unwrap()
}
