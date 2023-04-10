use rug::{Integer, integer::{Order, IsPrime}};
use rand::{CryptoRng, RngCore};

pub fn gen_prime<R: RngCore + CryptoRng>(rng: &mut R, num_bits: usize) -> Integer {
    loop {
        let p = sample_bits(rng, num_bits, true, true);
        if p.is_probably_prime(40) != IsPrime::No {
            println!("Found prime");
            return p;
        }
    }
}

pub fn gen_safe_prime<R: RngCore + CryptoRng>(rng: &mut R, num_bits: usize) -> Integer {
    loop {
        let p_half = gen_prime(rng, num_bits - 1);
        let p: Integer = 2 * p_half + 1;
        if p.is_probably_prime(40) != IsPrime::No {
            println!("Found safe prime");
            return p;
        }
    }
}

// Sample an integer from [0, 2**n). If exact_bits == true, sample from [2**(n-1), 2**n) instead. If
// odd is set, only sample odd integers.
fn sample_bits<R: RngCore + CryptoRng>(rng: &mut R, bits: usize, exact_bits: bool, odd: bool)
    -> Integer {
    let num_bytes = (bits + 7) / 8;
    let mut bytes = vec![0u8; num_bytes];
    rng.fill_bytes(&mut bytes);

    let high_bit = 1 << (((bits + 7) % 8) as u8);
    bytes[num_bytes - 1] &= high_bit + (high_bit - 1);
    if exact_bits {
        bytes[num_bytes - 1] |= high_bit;
    }
    if odd {
        bytes[0] |= 1;
    }

    Integer::from_digits(&bytes, Order::Lsf)
}
