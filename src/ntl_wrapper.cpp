#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include <memory>
#include <gmp.h>
#include <NTL/ZZ.h>

using namespace NTL;

static void ZZ_to_mpz(mpz_t output, const ZZ& a);

extern "C"
{
	void gen_prime_ntl(mpz_t output, const unsigned char* seed, size_t num_bits)
	{
		assert(num_bits <= LONG_MAX);
		assert(num_bits > 2);

		ZZ p;
		{
			RandomStreamPush push;
			SetSeed(seed, 32);
			GenPrime(p, num_bits);
		}

		return ZZ_to_mpz(output, p);
	}

	void gen_germain_prime_ntl(mpz_t output, const unsigned char* seed, size_t num_bits)
	{
		assert(num_bits <= LONG_MAX);
		assert(num_bits > 2);

		ZZ p;
		{
			RandomStreamPush push;
			SetSeed(seed, 32);
			GenGermainPrime(p, num_bits);
		}

		return ZZ_to_mpz(output, p);
	}
}

static void ZZ_to_mpz(mpz_t output, const ZZ& a)
{
	long num_bytes = NumBytes(a);
	std::unique_ptr<unsigned char[]> buf(new unsigned char[num_bytes]);
	BytesFromZZ(buf.get(), a, num_bytes);

	mpz_import(output, num_bytes, -1, 1, 0, CHAR_BIT - 8, buf.get());
}
