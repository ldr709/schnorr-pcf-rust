use core::mem::MaybeUninit;
use gmp_mpfr_sys::gmp::{bitcnt_t, mpz_t, mpz_init2};
use rug::Integer;
use rand::{CryptoRng, RngCore};

extern {
    fn gen_prime_ntl(output: *mut mpz_t, seed: *const u8, num_bits: usize);
    fn gen_germain_prime_ntl(output: *mut mpz_t, seed: *const u8, num_bits: usize);
}

pub fn gen_prime<R: RngCore + CryptoRng>(rng: &mut R, num_bits: usize) -> Integer {
    let mut seed = [0u8; 32];
    rng.fill_bytes(&mut seed);

    unsafe {
        let mut p_mpz = MaybeUninit::uninit();
        mpz_init2(p_mpz.as_mut_ptr(), num_bits as bitcnt_t);
        let mut p_mpz = p_mpz.assume_init();

        gen_prime_ntl(&mut p_mpz, seed.as_ptr(), num_bits);
        Integer::from_raw(p_mpz)
    }
}

pub fn gen_safe_prime<R: RngCore + CryptoRng>(rng: &mut R, num_bits: usize) -> Integer {
    let mut seed = [0u8; 32];
    rng.fill_bytes(&mut seed);

    let p_half = unsafe {
        let mut p_mpz = MaybeUninit::uninit();
        mpz_init2(p_mpz.as_mut_ptr(), num_bits as bitcnt_t);
        let mut p_mpz = p_mpz.assume_init();

        gen_germain_prime_ntl(&mut p_mpz, seed.as_ptr(), num_bits - 1);
        Integer::from_raw(p_mpz)
    };
    2 * p_half + 1
}


#[cfg(test)]
mod tests {
    use super::*;
    use rand_chacha::ChaCha20Rng;
    use rand::SeedableRng;
    use rug::integer::IsPrime;

    #[test]
    fn test_gen_prime() {
        let test_sizes = [
            (17, 10000),
            (23, 10000),
            (68, 5000),
            (129, 2000),
            (259, 500),
            (519, 100),
            (1536, 10),
            (2027, 5)
        ];

        let mut rng = ChaCha20Rng::from_entropy();
        for (len, cnt) in test_sizes {
            for _i in 0..cnt {
                let p = gen_prime(&mut rng, len);
                assert_eq!(p.significant_bits() as usize, len);
                assert_ne!(p.is_probably_prime(128), IsPrime::No);
            }
        }
    }

    #[test]
    fn test_gen_safe_prime() {
        let test_sizes = [
            (17, 1000),
            (23, 1000),
            (68, 500),
            (129, 200),
            (259, 50),
            (519, 10),
            (1536, 1),
        ];

        let mut rng = ChaCha20Rng::from_entropy();
        for (len, cnt) in test_sizes {
            for _i in 0..cnt {
                let p = gen_safe_prime(&mut rng, len);
                let p_half = (&p - Integer::from(1)) / Integer::from(2);
                assert_eq!(p.significant_bits() as usize, len);
                assert_ne!(p.is_probably_prime(128), IsPrime::No);
                assert_ne!(p_half.is_probably_prime(128), IsPrime::No);
            }
        }
    }
}
