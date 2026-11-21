#include "pch.h"
#include "convenience.h"

void sha::sha256_transform(uint32_t state[8], const uint8_t block[64])
{
	uint32_t w[64];
	uint32_t a, b, c, d, e, f, g, h;

	for (int i = 0; i < 16; i++)
	{
		w[i] = (block[i * 4] << 24) | (block[i * 4 + 1] << 16) |
			(block[i * 4 + 2] << 8) | (block[i * 4 + 3]);
	}

	for (int i = 16; i < 64; i++)
	{
		w[i] = SIG1(w[i - 2]) + w[i - 7] + SIG0(w[i - 15]) + w[i - 16];
	}

	a = state[0];
	b = state[1];
	c = state[2];
	d = state[3];
	e = state[4];
	f = state[5];
	g = state[6];
	h = state[7];

	for (int i = 0; i < 64; i++)
	{
		uint32_t temp1 = h + EP1(e) + CH(e, f, g) + k[i] + w[i];
		uint32_t temp2 = EP0(a) + MAJ(a, b, c);
		h = g;
		g = f;
		f = e;
		e = d + temp1;
		d = c;
		c = b;
		b = a;
		a = temp1 + temp2;
	}

	state[0] += a;
	state[1] += b;
	state[2] += c;
	state[3] += d;
	state[4] += e;
	state[5] += f;
	state[6] += g;
	state[7] += h;
}

// SHA-256 update function
void sha::sha256_update(uint32_t state[8], uint8_t buffer[64], const uint8_t *data, size_t len, uint64_t &bitlen)
{
	for (size_t i = 0; i < len; i++)
	{
		size_t buf_idx = (bitlen / 8) % 64;
		if (buf_idx >= 64) {
			std::cerr << "Buffer overflow detected in sha256_update!" << std::endl;
			return;
		}
		buffer[buf_idx] = data[i];
		bitlen += 8;
		if (bitlen % 512 == 0)
		{
			sha256_transform(state, buffer);
		}
	}
}

// SHA-256 padding and final hash computation
void sha::sha256_final(uint32_t state[8], uint8_t buffer[64], uint64_t bitlen, uint8_t hash[32])
{
	size_t i = bitlen / 8 % 64;

	if (i >= 64) {
		std::cerr << "Buffer overflow detected in sha256_final (initial)!" << std::endl;
		return;
	}
	buffer[i++] = 0x80;
	if (i > 56)
	{
		while (i < 64)
		{
			if (i >= 64) {
				std::cerr << "Buffer overflow detected in sha256_final (padding)!" << std::endl;
				return;
			}
			buffer[i++] = 0x00;
		}
		sha256_transform(state, buffer);
		i = 0;
	}

	while (i < 56)
	{
		if (i >= 64) {
			std::cerr << "Buffer overflow detected in sha256_final (final padding)!" << std::endl;
			return;
		}
		buffer[i++] = 0x00;
	}

	// memcpy to buffer + 56 is safe as buffer is 64 bytes
	bitlen = custom_bswap_64(bitlen);
	memcpy(buffer + 56, &bitlen, 8);
	sha256_transform(state, buffer);

	for (i = 0; i < 8; i++)
	{
		constexpr size_t HASH_SIZE = 32;
		const size_t off = static_cast<size_t>(i) * 4;
		if (off + 4 > HASH_SIZE) {
			std::cerr << "Buffer overflow detected in sha256_final (hash write)!" << std::endl;
			return;
		}
		// convert to big-endian value locally and copy 4 bytes
		uint32_t be = custom_bswap_32(state[i]);
		std::memcpy(hash + off, &be, sizeof(be));
	}
}

// Function to calculate SHA-256 hash
std::string sha::sha256(const std::string &input)
{
	uint32_t state[8] = {
		0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
		0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };

	uint8_t buffer[64] = { 0 };
	uint8_t hash[32];
	uint64_t bitlen = 0;

	sha256_update(state, buffer, reinterpret_cast<const uint8_t *>(input.c_str()), input.length(), bitlen);
	sha256_final(state, buffer, bitlen, hash);

	std::stringstream ss;
	for (int i = 0; i < 32; i++)
	{
		ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
	}

	return ss.str();
}
