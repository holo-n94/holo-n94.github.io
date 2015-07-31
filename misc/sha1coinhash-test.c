/*
 * gcc -O2 -o sha1coinhash-test sha1coinhash-test.c
 * gcc -O2 -msse2 -o sha1coinhash-test-sse2 sha1coinhash-test.c
 * gcc -O2 -msse2 -mxop -o sha1coinhash-test-xop sha1coinhash-test.c
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __XOP__
#include <x86intrin.h>
#endif


#if ((__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3))
#define WANT_BUILTIN_BSWAP
#else
#define bswap_32(x) ((((x) << 24) & 0xff000000u) | (((x) << 8) & 0x00ff0000u) \
                   | (((x) >> 8) & 0x0000ff00u) | (((x) >> 24) & 0x000000ffu))
#endif

static inline uint32_t swab32(uint32_t v)
{
#ifdef WANT_BUILTIN_BSWAP
	return __builtin_bswap32(v);
#else
	return bswap_32(v);
#endif
}


// constants and initial values defined in SHA-1
#define K0 0x5A827999
#define K1 0x6ED9EBA1
#define K2 0x8F1BBCDC
#define K3 0xCA62C1D6

#define H0 0x67452301
#define H1 0xEFCDAB89
#define H2 0x98BADCFE
#define H3 0x10325476
#define H4 0xC3D2E1F0


#define ROL32(_val32, _nBits) (((_val32)<<(_nBits))|((_val32)>>(32-(_nBits))))
#define Ch(x,y,z) ((x&(y^z))^z)
#define Maj(x,y,z) (((x|y)&z)|(x&y))

// W[t] = ROL32(W[t-3] ^ W[t-8] ^ W[t-14] ^ W[t-16], 1);
#define SHABLK(t) (W[t&15] = ROL32(W[(t+13)&15] ^ W[(t+8)&15] ^ W[(t+2)&15] ^ W[t&15], 1))

#define _RS0(v,w,x,y,z,i) { z += Ch(w,x,y) + i + K0 + ROL32(v,5);  w=ROL32(w,30); }
#define _RS00(v,w,x,y,z)  { z += Ch(w,x,y) + K0 + ROL32(v,5);  w=ROL32(w,30); }
#define _RS1(v,w,x,y,z,i) { z += (w^x^y) + i + K1 + ROL32(v,5);  w=ROL32(w,30); }

#define _R0(v,w,x,y,z,t) { z += Ch(w,x,y) + SHABLK(t) + K0 + ROL32(v,5);  w=ROL32(w,30); }
#define _R1(v,w,x,y,z,t) { z += (w^x^y) + SHABLK(t) + K1 + ROL32(v,5);  w=ROL32(w,30); }
#define _R2(v,w,x,y,z,t) { z += Maj(w,x,y) + SHABLK(t) + K2 + ROL32(v,5);  w=ROL32(w,30); }
#define _R3(v,w,x,y,z,t) { z += (w^x^y) + SHABLK(t) + K3 + ROL32(v,5);  w=ROL32(w,30); }


void sha1hash12byte(const char *input, uint32_t *m_state)
{
	uint32_t W[16];
	uint32_t a, b, c, d, e;
	int i;

	// SHA-1 initialization constants
	m_state[0] = H0;
	m_state[1] = H1;
	m_state[2] = H2;
	m_state[3] = H3;
	m_state[4] = H4;

	a = m_state[0];
	b = m_state[1];
	c = m_state[2];
	d = m_state[3];
	e = m_state[4];

	// input[0] to input[11], 12byte, 96bits
	for (i = 0; i < 3; i++){
		W[i] = input[4*i+0] << 24 | input[4*i+1] << 16 | input[4*i+2] << 8 | input[4*i+3];
	}

	W[3] = 0x80000000;		// padding

/*
	for (int i = 4; i < 15; i++){
		W[i] = 0;
	}
*/

	W[15] = 12 * 8;		// bits of Message Block (12 bytes * 8 bits)

	// round 0 to 15
	_RS0(a, b, c, d, e, W[0]);
	_RS0(e, a, b, c, d, W[1]);
	_RS0(d, e, a, b, c, W[2]);
	_RS0(c, d, e, a, b, W[3]);
	_RS00(b, c, d, e, a);		// W[4] == 0
	_RS00(a, b, c, d, e);		// W[5] == 0
	_RS00(e, a, b, c, d);		// W[6] == 0
	_RS00(d, e, a, b, c);		// W[7] == 0
	_RS00(c, d, e, a, b);		// W[8] == 0
	_RS00(b, c, d, e, a);		// W[9] == 0
	_RS00(a, b, c, d, e);		// W[10] == 0
	_RS00(e, a, b, c, d);		// W[11] == 0
	_RS00(d, e, a, b, c);		// W[12] == 0
	_RS00(c, d, e, a, b);		// W[13] == 0
	_RS00(b, c, d, e, a);		// W[14] == 0
	_RS0(a, b, c, d, e, W[15]);

	// round 16 to 19
	W[0] = ROL32(W[2] ^ W[0], 1);		// (t, W[t-3], W[t-8], W[t-14], W[t-16]) = (16, W[13]==0, W[8]==0, W[2], W[0])
	_RS0(e, a, b, c, d, W[0]);

	W[1] = ROL32(W[3] ^ W[1], 1);		// (17, W[14]==0, W[9]==0, W[3], W[1])
	_RS0(d, e, a, b, c, W[1]);

	W[2] = ROL32(W[15] ^ W[2], 1);		// (18, W[15], W[10]==0, W[4]==0, W[2])
	_RS0(c, d, e, a, b, W[2]);

	W[3] = ROL32(W[0] ^ W[3], 1);		// (19, W[0], W[11]==0, W[5]==0, W[3])
	_RS0(b, c, d, e, a, W[3]);

	// round 20 to 31
	W[4] = ROL32(W[1], 1);				// (20, W[1], W[12]==0, W[6]==0, W[4]==0)
	_RS1(a, b, c, d, e, W[4]);

	W[5] = ROL32(W[2], 1);				// (21, W[2], W[13]==0, W[7]==0, W[5]==0)
	_RS1(e, a, b, c, d, W[5]);

	W[6] = ROL32(W[3], 1);				// (22, W[3], W[14]==0, W[8]==0, W[6]==0)
	_RS1(d, e, a, b, c, W[6]);

	W[7] = ROL32(W[4] ^ W[15], 1);		// (23, W[4], W[15], W[9]==0, W[7]==0)
	_RS1(c, d, e, a, b, W[7]);

	W[8] = ROL32(W[5] ^ W[0], 1);		// (24, W[5], W[0], W[10]==0, W[8]==0)
	_RS1(b, c, d, e, a, W[8]);

	W[9] = ROL32(W[6] ^ W[1], 1);		// (25, W[6], W[1], W[11]==0, W[9]==0)
	_RS1(a, b, c, d, e, W[9]);

	W[10] = ROL32(W[7] ^ W[2], 1);		// (26, W[7], W[2], W[12]==0, W[10]==0)
	_RS1(e, a, b, c, d, W[10]);

	W[11] = ROL32(W[8] ^ W[3], 1);		// (27, W[8], W[3], W[13]==0, W[11]==0)
	_RS1(d, e, a, b, c, W[11]);

	W[12] = ROL32(W[9] ^ W[4], 1);		// (28, W[9], W[4], W[14]==0, W[12]==0)
	_RS1(c, d, e, a, b, W[12]);

	W[13] = ROL32(W[10] ^ W[5] ^ W[15], 1);		// (29, W[10], W[5], W[15], W[13]==0)
	_RS1(b, c, d, e, a, W[13]);

	W[14] = ROL32(W[11] ^ W[6] ^ W[0], 1);		// (30, W[11], W[6], W[0], W[14]==0)
	_RS1(a, b, c, d, e, W[14]);

	W[15] = ROL32(W[12] ^ W[7] ^ W[1] ^ W[15], 1);		// (31, W[12], W[7], W[1], W[15])
	_RS1(e, a, b, c, d, W[15]);

	// round 32 to 39
	_R1(d, e, a, b, c, 32);
	_R1(c, d, e, a, b, 33);
	_R1(b, c, d, e, a, 34);
	_R1(a, b, c, d, e, 35);
	_R1(e, a, b, c, d, 36);
	_R1(d, e, a, b, c, 37);
	_R1(c, d, e, a, b, 38);
	_R1(b, c, d, e, a, 39);

	// round 40 to 59
	_R2(a, b, c, d, e, 40);
	_R2(e, a, b, c, d, 41);
	_R2(d, e, a, b, c, 42);
	_R2(c, d, e, a, b, 43);
	_R2(b, c, d, e, a, 44);
	_R2(a, b, c, d, e, 45);
	_R2(e, a, b, c, d, 46);
	_R2(d, e, a, b, c, 47);
	_R2(c, d, e, a, b, 48);
	_R2(b, c, d, e, a, 49);
	_R2(a, b, c, d, e, 50);
	_R2(e, a, b, c, d, 51);
	_R2(d, e, a, b, c, 52);
	_R2(c, d, e, a, b, 53);
	_R2(b, c, d, e, a, 54);
	_R2(a, b, c, d, e, 55);
	_R2(e, a, b, c, d, 56);
	_R2(d, e, a, b, c, 57);
	_R2(c, d, e, a, b, 58);
	_R2(b, c, d, e, a, 59);

	// round 60 to 79
	_R3(a, b, c, d, e, 60);
	_R3(e, a, b, c, d, 61);
	_R3(d, e, a, b, c, 62);
	_R3(c, d, e, a, b, 63);
	_R3(b, c, d, e, a, 64);
	_R3(a, b, c, d, e, 65);
	_R3(e, a, b, c, d, 66);
	_R3(d, e, a, b, c, 67);
	_R3(c, d, e, a, b, 68);
	_R3(b, c, d, e, a, 69);
	_R3(a, b, c, d, e, 70);
	_R3(e, a, b, c, d, 71);
	_R3(d, e, a, b, c, 72);
	_R3(c, d, e, a, b, 73);
	_R3(b, c, d, e, a, 74);
	_R3(a, b, c, d, e, 75);
	_R3(e, a, b, c, d, 76);
	_R3(d, e, a, b, c, 77);
	_R3(c, d, e, a, b, 78);
	_R3(b, c, d, e, a, 79);

	// Add the working vars back into state
	m_state[0] += a;
	m_state[1] += b;
	m_state[2] += c;
	m_state[3] += d;
	m_state[4] += e;

	m_state[0] = swab32(m_state[0]);
	m_state[1] = swab32(m_state[1]);
	m_state[2] = swab32(m_state[2]);
	m_state[3] = swab32(m_state[3]);
	m_state[4] = swab32(m_state[4]);
}


void sha1hash80byte(const uint8_t *input, uint32_t *m_state)
{
	uint32_t W[16];
	uint32_t a, b, c, d, e;
	int i;

	// SHA-1 initialization constants
	m_state[0] = H0;
	m_state[1] = H1;
	m_state[2] = H2;
	m_state[3] = H3;
	m_state[4] = H4;

	a = m_state[0];
	b = m_state[1];
	c = m_state[2];
	d = m_state[3];
	e = m_state[4];

	// input[0] to input[63], 64bytes, 512bits
	for (i = 0; i < 16; i++){
		W[i] = input[4*i+0] << 24 | input[4*i+1] << 16 | input[4*i+2] << 8 | input[4*i+3];
	}

	// round 0 to 15
	_RS0(a, b, c, d, e, W[0]);
	_RS0(e, a, b, c, d, W[1]);
	_RS0(d, e, a, b, c, W[2]);
	_RS0(c, d, e, a, b, W[3]);
	_RS0(b, c, d, e, a, W[4]);
	_RS0(a, b, c, d, e, W[5]);
	_RS0(e, a, b, c, d, W[6]);
	_RS0(d, e, a, b, c, W[7]);
	_RS0(c, d, e, a, b, W[8]);
	_RS0(b, c, d, e, a, W[9]);
	_RS0(a, b, c, d, e, W[10]);
	_RS0(e, a, b, c, d, W[11]);
	_RS0(d, e, a, b, c, W[12]);
	_RS0(c, d, e, a, b, W[13]);
	_RS0(b, c, d, e, a, W[14]);
	_RS0(a, b, c, d, e, W[15]);

	// round 16 to 19
	_R0(e, a, b, c, d, 16);
	_R0(d, e, a, b, c, 17);
	_R0(c, d, e, a, b, 18);
	_R0(b, c, d, e, a, 19);

	// round 20 to 39
	_R1(a, b, c, d, e, 20);
	_R1(e, a, b, c, d, 21);
	_R1(d, e, a, b, c, 22);
	_R1(c, d, e, a, b, 23);
	_R1(b, c, d, e, a, 24);
	_R1(a, b, c, d, e, 25);
	_R1(e, a, b, c, d, 26);
	_R1(d, e, a, b, c, 27);
	_R1(c, d, e, a, b, 28);
	_R1(b, c, d, e, a, 29);
	_R1(a, b, c, d, e, 30);
	_R1(e, a, b, c, d, 31);
	_R1(d, e, a, b, c, 32);
	_R1(c, d, e, a, b, 33);
	_R1(b, c, d, e, a, 34);
	_R1(a, b, c, d, e, 35);
	_R1(e, a, b, c, d, 36);
	_R1(d, e, a, b, c, 37);
	_R1(c, d, e, a, b, 38);
	_R1(b, c, d, e, a, 39);

	// round 40 to 59
	_R2(a, b, c, d, e, 40);
	_R2(e, a, b, c, d, 41);
	_R2(d, e, a, b, c, 42);
	_R2(c, d, e, a, b, 43);
	_R2(b, c, d, e, a, 44);
	_R2(a, b, c, d, e, 45);
	_R2(e, a, b, c, d, 46);
	_R2(d, e, a, b, c, 47);
	_R2(c, d, e, a, b, 48);
	_R2(b, c, d, e, a, 49);
	_R2(a, b, c, d, e, 50);
	_R2(e, a, b, c, d, 51);
	_R2(d, e, a, b, c, 52);
	_R2(c, d, e, a, b, 53);
	_R2(b, c, d, e, a, 54);
	_R2(a, b, c, d, e, 55);
	_R2(e, a, b, c, d, 56);
	_R2(d, e, a, b, c, 57);
	_R2(c, d, e, a, b, 58);
	_R2(b, c, d, e, a, 59);

	// round 60 to 79
	_R3(a, b, c, d, e, 60);
	_R3(e, a, b, c, d, 61);
	_R3(d, e, a, b, c, 62);
	_R3(c, d, e, a, b, 63);
	_R3(b, c, d, e, a, 64);
	_R3(a, b, c, d, e, 65);
	_R3(e, a, b, c, d, 66);
	_R3(d, e, a, b, c, 67);
	_R3(c, d, e, a, b, 68);
	_R3(b, c, d, e, a, 69);
	_R3(a, b, c, d, e, 70);
	_R3(e, a, b, c, d, 71);
	_R3(d, e, a, b, c, 72);
	_R3(c, d, e, a, b, 73);
	_R3(b, c, d, e, a, 74);
	_R3(a, b, c, d, e, 75);
	_R3(e, a, b, c, d, 76);
	_R3(d, e, a, b, c, 77);
	_R3(c, d, e, a, b, 78);
	_R3(b, c, d, e, a, 79);

	// Add the working vars back into state
	m_state[0] += a;
	m_state[1] += b;
	m_state[2] += c;
	m_state[3] += d;
	m_state[4] += e;

	a = m_state[0];
	b = m_state[1];
	c = m_state[2];
	d = m_state[3];
	e = m_state[4];

	// input[64] to input[79], 16bytes, 128bits
	for (i = 0; i < 4; i++){
		W[i] = input[4*i+64] << 24 | input[4*i+65] << 16 | input[4*i+66] << 8 | input[4*i+67];
	}

	W[4] = 0x80000000;		// padding

/*
	for (int i = 5; i < 15; i++){
		W[i] = 0;
	}
*/

	W[15] = 80 * 8;		// bits of Message Block (80 bytes * 8 bits)

	// round 0 to 15
	_RS0(a, b, c, d, e, W[0]);
	_RS0(e, a, b, c, d, W[1]);
	_RS0(d, e, a, b, c, W[2]);
	_RS0(c, d, e, a, b, W[3]);
	_RS0(b, c, d, e, a, W[4]);
	_RS00(a, b, c, d, e);		// W[5] == 0
	_RS00(e, a, b, c, d);		// W[6] == 0
	_RS00(d, e, a, b, c);		// W[7] == 0
	_RS00(c, d, e, a, b);		// W[8] == 0
	_RS00(b, c, d, e, a);		// W[9] == 0
	_RS00(a, b, c, d, e);		// W[10] == 0
	_RS00(e, a, b, c, d);		// W[11] == 0
	_RS00(d, e, a, b, c);		// W[12] == 0
	_RS00(c, d, e, a, b);		// W[13] == 0
	_RS00(b, c, d, e, a);		// W[14] == 0
	_RS0(a, b, c, d, e, W[15]);

	// round 16 to 19
	W[0] = ROL32(W[2] ^ W[0], 1);		// (t, W[t-3], W[t-8], W[t-14], W[t-16]) = (16, W[13]==0, W[8]==0, W[2], W[0])
	_RS0(e, a, b, c, d, W[0]);

	W[1] = ROL32(W[3] ^ W[1], 1);		// (17, W[14]==0, W[9]==0, W[3], W[1])
	_RS0(d, e, a, b, c, W[1]);

	W[2] = ROL32(W[15] ^ W[4] ^ W[2], 1);		// (18, W[15], W[10]==0, W[4], W[2])
	_RS0(c, d, e, a, b, W[2]);

	W[3] = ROL32(W[0] ^ W[3], 1);		// (19, W[0], W[11]==0, W[5]==0, W[3])
	_RS0(b, c, d, e, a, W[3]);

	// round 20 to 31
	W[4] = ROL32(W[1] ^ W[4], 1);		// (20, W[1], W[12]==0, W[6]==0, W[4])
	_RS1(a, b, c, d, e, W[4]);

	W[5] = ROL32(W[2], 1);				// (21, W[2], W[13]==0, W[7]==0, W[5]==0)
	_RS1(e, a, b, c, d, W[5]);

	W[6] = ROL32(W[3], 1);				// (22, W[3], W[14]==0, W[8]==0, W[6]==0)
	_RS1(d, e, a, b, c, W[6]);

	W[7] = ROL32(W[4] ^ W[15], 1);		// (23, W[4], W[15], W[9]==0, W[7]==0)
	_RS1(c, d, e, a, b, W[7]);

	W[8] = ROL32(W[5] ^ W[0], 1);		// (24, W[5], W[0], W[10]==0, W[8]==0)
	_RS1(b, c, d, e, a, W[8]);

	W[9] = ROL32(W[6] ^ W[1], 1);		// (25, W[6], W[1], W[11]==0, W[9]==0)
	_RS1(a, b, c, d, e, W[9]);

	W[10] = ROL32(W[7] ^ W[2], 1);		// (26, W[7], W[2], W[12]==0, W[10]==0)
	_RS1(e, a, b, c, d, W[10]);

	W[11] = ROL32(W[8] ^ W[3], 1);		// (27, W[8], W[3], W[13]==0, W[11]==0)
	_RS1(d, e, a, b, c, W[11]);

	W[12] = ROL32(W[9] ^ W[4], 1);		// (28, W[9], W[4], W[14]==0, W[12]==0)
	_RS1(c, d, e, a, b, W[12]);

	W[13] = ROL32(W[10] ^ W[5] ^ W[15], 1);		// (29, W[10], W[5], W[15], W[13]==0)
	_RS1(b, c, d, e, a, W[13]);

	W[14] = ROL32(W[11] ^ W[6] ^ W[0], 1);		// (30, W[11], W[6], W[0], W[14]==0)
	_RS1(a, b, c, d, e, W[14]);

	W[15] = ROL32(W[12] ^ W[7] ^ W[1] ^ W[15], 1);		// (31, W[12], W[7], W[1], W[15])
	_RS1(e, a, b, c, d, W[15]);

	// round 32 to 39
	_R1(d, e, a, b, c, 32);
	_R1(c, d, e, a, b, 33);
	_R1(b, c, d, e, a, 34);
	_R1(a, b, c, d, e, 35);
	_R1(e, a, b, c, d, 36);
	_R1(d, e, a, b, c, 37);
	_R1(c, d, e, a, b, 38);
	_R1(b, c, d, e, a, 39);

	// round 40 to 59
	_R2(a, b, c, d, e, 40);
	_R2(e, a, b, c, d, 41);
	_R2(d, e, a, b, c, 42);
	_R2(c, d, e, a, b, 43);
	_R2(b, c, d, e, a, 44);
	_R2(a, b, c, d, e, 45);
	_R2(e, a, b, c, d, 46);
	_R2(d, e, a, b, c, 47);
	_R2(c, d, e, a, b, 48);
	_R2(b, c, d, e, a, 49);
	_R2(a, b, c, d, e, 50);
	_R2(e, a, b, c, d, 51);
	_R2(d, e, a, b, c, 52);
	_R2(c, d, e, a, b, 53);
	_R2(b, c, d, e, a, 54);
	_R2(a, b, c, d, e, 55);
	_R2(e, a, b, c, d, 56);
	_R2(d, e, a, b, c, 57);
	_R2(c, d, e, a, b, 58);
	_R2(b, c, d, e, a, 59);

	// round 60 to 79
	_R3(a, b, c, d, e, 60);
	_R3(e, a, b, c, d, 61);
	_R3(d, e, a, b, c, 62);
	_R3(c, d, e, a, b, 63);
	_R3(b, c, d, e, a, 64);
	_R3(a, b, c, d, e, 65);
	_R3(e, a, b, c, d, 66);
	_R3(d, e, a, b, c, 67);
	_R3(c, d, e, a, b, 68);
	_R3(b, c, d, e, a, 69);
	_R3(a, b, c, d, e, 70);
	_R3(e, a, b, c, d, 71);
	_R3(d, e, a, b, c, 72);
	_R3(c, d, e, a, b, 73);
	_R3(b, c, d, e, a, 74);
	_R3(a, b, c, d, e, 75);
	_R3(e, a, b, c, d, 76);
	_R3(d, e, a, b, c, 77);
	_R3(c, d, e, a, b, 78);
	_R3(b, c, d, e, a, 79);

	// Add the working vars back into state
	m_state[0] += a;
	m_state[1] += b;
	m_state[2] += c;
	m_state[3] += d;
	m_state[4] += e;

	m_state[0] = swab32(m_state[0]);
	m_state[1] = swab32(m_state[1]);
	m_state[2] = swab32(m_state[2]);
	m_state[3] = swab32(m_state[3]);
	m_state[4] = swab32(m_state[4]);
}


void b64enc(const uint32_t *hash, char *str)
{
	const char b64t[] = {
		'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
		'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
		'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
		'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '/'
	};

	str[0] = b64t[hash[0] >> 26];
	str[1] = b64t[(hash[0] >> 20) & 63];
	str[2] = b64t[(hash[0] >> 14) & 63];
	str[3] = b64t[(hash[0] >> 8) & 63];
	str[4] = b64t[(hash[0] >> 2) & 63];
	str[5] = b64t[(hash[0] << 4 | hash[1] >> 28) & 63];
	str[6] = b64t[(hash[1] >> 22) & 63];
	str[7] = b64t[(hash[1] >> 16) & 63];
	str[8] = b64t[(hash[1] >> 10) & 63];
	str[9] = b64t[(hash[1] >> 4) & 63];
	str[10] = b64t[(hash[1] << 2 | hash[2] >> 30) & 63];
	str[11] = b64t[(hash[2] >> 24) & 63];
	str[12] = b64t[(hash[2] >> 18) & 63];
	str[13] = b64t[(hash[2] >> 12) & 63];
	str[14] = b64t[(hash[2] >> 6) & 63];
	str[15] = b64t[hash[2] & 63];
	str[16] = b64t[hash[3] >> 26];
	str[17] = b64t[(hash[3] >> 20) & 63];
	str[18] = b64t[(hash[3] >> 14) & 63];
	str[19] = b64t[(hash[3] >> 8) & 63];
	str[20] = b64t[(hash[3] >> 2) & 63];
	str[21] = b64t[(hash[3] << 4 | hash[4] >> 28) & 63];
	str[22] = b64t[(hash[4] >> 22) & 63];
	str[23] = b64t[(hash[4] >> 16) & 63];
	str[24] = b64t[(hash[4] >> 10) & 63];
	str[25] = b64t[(hash[4] >> 4) & 63];
	str[26] = b64t[(hash[4] << 2) & 63];
	str[27] = 0;
}


#ifdef __SSE2__		//////////////////////////////////////////////////

#define MM_OR(a, b) _mm_or_si128((a), (b))
#define MM_AND(a, b) _mm_and_si128((a), (b))
#define MM_XOR(a, b) _mm_xor_si128((a), (b))
#define MM_ADD(a, b) _mm_add_epi32((a), (b))

#define MM_SLLI(a, b) _mm_slli_epi32((a), (b))
#define MM_SRLI(a, b) _mm_srli_epi32((a), (b))

#define MM_SET1(a) _mm_set1_epi32((a))

#define MM_LOAD(a) _mm_load_si128((__m128i *)(a))
#define MM_STORE(a, b) _mm_store_si128((__m128i *)(a), (b))

#undef ROL32
#undef Ch
#undef Maj

#ifdef __XOP__
#define ROL32(_val32, _nBits) _mm_roti_epi32((_val32), (_nBits))
#define Ch(x,y,z) _mm_cmov_si128((y),(z),(x))
#define Maj(x,y,z) Ch(MM_XOR((x),(z)),(y),(z))
#else
#define ROL32(_val32, _nBits) (MM_OR(MM_SLLI((_val32), (_nBits)), MM_SRLI((_val32), 32-(_nBits))))
#define Ch(x,y,z) (MM_XOR(MM_AND((x), MM_XOR((y), (z))), (z)))
#define Maj(x,y,z) (MM_OR(MM_AND(MM_OR((x), (y)), (z)), MM_AND((x), (y))))
#endif

#undef SHABLK
#define SHABLK(t) (W[(t)&15] = ROL32(MM_XOR(MM_XOR(MM_XOR(W[((t)+13)&15], W[((t)+8)&15]), W[((t)+2)&15]), W[(t)&15]), 1))

#undef _RS0
#define _RS0(v,w,x,y,z,i) { \
	z = MM_ADD((z), MM_ADD(MM_ADD(MM_ADD(Ch(w,x,y), (i)), MM_SET1(K0)), ROL32(v,5))); \
	w = ROL32(w, 30); \
}

#undef _RS00
#define _RS00(v,w,x,y,z) { \
	z = MM_ADD((z), MM_ADD(MM_ADD(Ch(w,x,y), MM_SET1(K0)), ROL32(v,5))); \
	w = ROL32(w, 30); \
}

#undef _RS1
#define _RS1(v,w,x,y,z,i) { \
	z = MM_ADD((z), MM_ADD(MM_ADD(MM_ADD(MM_XOR(MM_XOR((w), (x)), (y)), (i)), MM_SET1(K1)), ROL32(v,5))); \
	w = ROL32(w, 30); \
}

#undef _R0
#define _R0(v,w,x,y,z,t) { \
	z = MM_ADD((z), MM_ADD(MM_ADD(MM_ADD(Ch(w,x,y), SHABLK(t)), MM_SET1(K0)), ROL32(v,5))); \
	w = ROL32(w, 30); \
}

#undef _R1
#define _R1(v,w,x,y,z,t) { \
	z = MM_ADD((z), MM_ADD(MM_ADD(MM_ADD(MM_XOR(MM_XOR((w), (x)), (y)), SHABLK(t)), MM_SET1(K1)), ROL32(v,5))); \
	w = ROL32(w, 30); \
}

#undef _R2
#define _R2(v,w,x,y,z,t) { \
	z = MM_ADD((z), MM_ADD(MM_ADD(MM_ADD(Maj(w,x,y), SHABLK(t)), MM_SET1(K2)), ROL32(v,5))); \
	w = ROL32(w, 30); \
}

#undef _R3
#define _R3(v,w,x,y,z,t) { \
	z = MM_ADD((z), MM_ADD(MM_ADD(MM_ADD(MM_XOR(MM_XOR((w), (x)), (y)), SHABLK(t)), MM_SET1(K3)), ROL32(v,5))); \
	w = ROL32(w, 30); \
}


void sha1hash12byte_sse2(const char *input, __m128i *m_state)
{
	__attribute__((aligned(16))) __m128i W[16];
	__attribute__((aligned(16))) __m128i a, b, c, d, e;
	int i;

	// SHA-1 initialization constants
	m_state[0] = MM_SET1(H0);
	m_state[1] = MM_SET1(H1);
	m_state[2] = MM_SET1(H2);
	m_state[3] = MM_SET1(H3);
	m_state[4] = MM_SET1(H4);

	a = m_state[0];
	b = m_state[1];
	c = m_state[2];
	d = m_state[3];
	e = m_state[4];

	for (i = 0; i < 3; i++){
		W[i] = MM_LOAD(&input[16 * i]);
	}

	W[3] = MM_SET1(0x80000000);		// padding

	W[15] = MM_SET1(12 * 8);		// bits of Message Block (12 bytes * 8 bits)

	// round 0 to 15
	_RS0(a, b, c, d, e, W[0]);
	_RS0(e, a, b, c, d, W[1]);
	_RS0(d, e, a, b, c, W[2]);
	_RS0(c, d, e, a, b, W[3]);
	_RS00(b, c, d, e, a);		// W[4] == 0
	_RS00(a, b, c, d, e);		// W[5] == 0
	_RS00(e, a, b, c, d);		// W[6] == 0
	_RS00(d, e, a, b, c);		// W[7] == 0
	_RS00(c, d, e, a, b);		// W[8] == 0
	_RS00(b, c, d, e, a);		// W[9] == 0
	_RS00(a, b, c, d, e);		// W[10] == 0
	_RS00(e, a, b, c, d);		// W[11] == 0
	_RS00(d, e, a, b, c);		// W[12] == 0
	_RS00(c, d, e, a, b);		// W[13] == 0
	_RS00(b, c, d, e, a);		// W[14] == 0
	_RS0(a, b, c, d, e, W[15]);

	// round 16 to 19
	// (t, W[t-3], W[t-8], W[t-14], W[t-16]) = (16, W[13]==0, W[8]==0, W[2], W[0])
	W[0] = ROL32(MM_XOR(W[2], W[0]), 1);
	_RS0(e, a, b, c, d, W[0]);

	// (17, W[14]==0, W[9]==0, W[3], W[1])
	W[1] = ROL32(MM_XOR(W[3], W[1]), 1);
	_RS0(d, e, a, b, c, W[1]);

	// (18, W[15], W[10]==0, W[4]==0, W[2])
	W[2] = ROL32(MM_XOR(W[15], W[2]), 1);
	_RS0(c, d, e, a, b, W[2]);

	// (19, W[0], W[11]==0, W[5]==0, W[3])
	W[3] = ROL32(MM_XOR(W[0], W[3]), 1);
	_RS0(b, c, d, e, a, W[3]);

	// round 20 to 31
	// (20, W[1], W[12]==0, W[6]==0, W[4]==0)
	W[4] = ROL32(W[1], 1);
	_RS1(a, b, c, d, e, W[4]);

	// (21, W[2], W[13]==0, W[7]==0, W[5]==0)
	W[5] = ROL32(W[2], 1);
	_RS1(e, a, b, c, d, W[5]);

	// (22, W[3], W[14]==0, W[8]==0, W[6]==0)
	W[6] = ROL32(W[3], 1);
	_RS1(d, e, a, b, c, W[6]);

	// (23, W[4], W[15], W[9]==0, W[7]==0)
	W[7] = ROL32(MM_XOR(W[4], W[15]), 1);
	_RS1(c, d, e, a, b, W[7]);

	// (24, W[5], W[0], W[10]==0, W[8]==0)
	W[8] = ROL32(MM_XOR(W[5], W[0]), 1);
	_RS1(b, c, d, e, a, W[8]);

	// (25, W[6], W[1], W[11]==0, W[9]==0)
	W[9] = ROL32(MM_XOR(W[6], W[1]), 1);
	_RS1(a, b, c, d, e, W[9]);

	// (26, W[7], W[2], W[12]==0, W[10]==0)
	W[10] = ROL32(MM_XOR(W[7], W[2]), 1);
	_RS1(e, a, b, c, d, W[10]);

	// (27, W[8], W[3], W[13]==0, W[11]==0)
	W[11] = ROL32(MM_XOR(W[8], W[3]), 1);
	_RS1(d, e, a, b, c, W[11]);

	// (28, W[9], W[4], W[14]==0, W[12]==0)
	W[12] = ROL32(MM_XOR(W[9], W[4]), 1);
	_RS1(c, d, e, a, b, W[12]);

	// (29, W[10], W[5], W[15], W[13]==0)
	W[13] = ROL32(MM_XOR(MM_XOR(W[10], W[5]), W[15]), 1);
	_RS1(b, c, d, e, a, W[13]);

	// (30, W[11], W[6], W[0], W[14]==0)
	W[14] = ROL32(MM_XOR(MM_XOR(W[11], W[6]), W[0]), 1);
	_RS1(a, b, c, d, e, W[14]);

	// (31, W[12], W[7], W[1], W[15])
	W[15] = ROL32(MM_XOR(MM_XOR(MM_XOR(W[12], W[7]), W[1]), W[15]), 1);
	_RS1(e, a, b, c, d, W[15]);

	// round 32 to 39
	_R1(d, e, a, b, c, 32);
	_R1(c, d, e, a, b, 33);
	_R1(b, c, d, e, a, 34);
	_R1(a, b, c, d, e, 35);
	_R1(e, a, b, c, d, 36);
	_R1(d, e, a, b, c, 37);
	_R1(c, d, e, a, b, 38);
	_R1(b, c, d, e, a, 39);

	// round 40 to 59
	_R2(a, b, c, d, e, 40);
	_R2(e, a, b, c, d, 41);
	_R2(d, e, a, b, c, 42);
	_R2(c, d, e, a, b, 43);
	_R2(b, c, d, e, a, 44);
	_R2(a, b, c, d, e, 45);
	_R2(e, a, b, c, d, 46);
	_R2(d, e, a, b, c, 47);
	_R2(c, d, e, a, b, 48);
	_R2(b, c, d, e, a, 49);
	_R2(a, b, c, d, e, 50);
	_R2(e, a, b, c, d, 51);
	_R2(d, e, a, b, c, 52);
	_R2(c, d, e, a, b, 53);
	_R2(b, c, d, e, a, 54);
	_R2(a, b, c, d, e, 55);
	_R2(e, a, b, c, d, 56);
	_R2(d, e, a, b, c, 57);
	_R2(c, d, e, a, b, 58);
	_R2(b, c, d, e, a, 59);

	// round 60 to 79
	_R3(a, b, c, d, e, 60);
	_R3(e, a, b, c, d, 61);
	_R3(d, e, a, b, c, 62);
	_R3(c, d, e, a, b, 63);
	_R3(b, c, d, e, a, 64);
	_R3(a, b, c, d, e, 65);
	_R3(e, a, b, c, d, 66);
	_R3(d, e, a, b, c, 67);
	_R3(c, d, e, a, b, 68);
	_R3(b, c, d, e, a, 69);
	_R3(a, b, c, d, e, 70);
	_R3(e, a, b, c, d, 71);
	_R3(d, e, a, b, c, 72);
	_R3(c, d, e, a, b, 73);
	_R3(b, c, d, e, a, 74);
	_R3(a, b, c, d, e, 75);
	_R3(e, a, b, c, d, 76);
	_R3(d, e, a, b, c, 77);
	_R3(c, d, e, a, b, 78);
	_R3(b, c, d, e, a, 79);

	// Add the working vars back into state
	m_state[0] = MM_ADD(m_state[0], a);
	m_state[1] = MM_ADD(m_state[1], b);
	m_state[2] = MM_ADD(m_state[2], c);
	m_state[3] = MM_ADD(m_state[3], d);
	m_state[4] = MM_ADD(m_state[4], e);
}


#endif	// __SSE2__		//////////////////////////////////////////////////


#ifdef __SSE2__		//////////////////////////////////////////////////
void sha1coinhash(uint8_t *input, uint32_t *hash)
{
	char str[40] __attribute__((aligned(32))) = {0}; // 26 + 11 + 1 + padding
	char tripkey[4 * 4 * 3] __attribute__((aligned(32)));

	uint32_t prehash[5] __attribute__((aligned(32)));
	__m128i prehash_m128i[5] __attribute__((aligned(32)));
	__m128i hash_m128i[5] __attribute__((aligned(32)));

	int i, j, k;
	__attribute__((aligned(16))) uint32_t tmp[4];

	sha1hash80byte(input, prehash);
	b64enc(prehash, str);
	memcpy(str + 26, str, 11);

	for (i = 0; i < 5; i++){
		hash_m128i[i] = _mm_set1_epi32(0);
	}

	for (k = 0; k < 7; k++){
		// generate tripkey table from str
		for (i = 0; i < 4; i++){
			for (j = 0; j < 3; j++){
				tripkey[3 + 4 * i + 16 * j] = str[0 + 4 * k + i + 4 * j];
				tripkey[2 + 4 * i + 16 * j] = str[1 + 4 * k + i + 4 * j];
				tripkey[1 + 4 * i + 16 * j] = str[2 + 4 * k + i + 4 * j];
				tripkey[0 + 4 * i + 16 * j] = str[3 + 4 * k + i + 4 * j];
			}
		}

		sha1hash12byte_sse2(tripkey, prehash_m128i);

		if (k < 6){
			for (i = 0; i < 5; i++){
				hash_m128i[i] = _mm_xor_si128(hash_m128i[i], prehash_m128i[i]);
			}
		}
	}

	for (i = 0; i < 5; i++){
		_mm_store_si128((__m128i *)tmp, hash_m128i[i]);
		hash[i] = tmp[0] ^ tmp[1] ^ tmp[2] ^ tmp[3];

		_mm_store_si128((__m128i *)tmp, prehash_m128i[i]);
		hash[i] ^= tmp[0] ^ tmp[1];

		hash[i] = swab32(hash[i]);
	}
}

#else	// ifudef __SSE2__		//////////////////////////////////////////////////

void sha1coinhash(uint8_t *input, uint32_t *hash)
{
	char str[38] __attribute__((aligned(32))); // 26 + 11 + 1
	uint32_t prehash[5] __attribute__((aligned(32)));
	int i;

	sha1hash80byte(input, prehash);
	b64enc(prehash, str);
	memcpy(str + 26, str, 11);

	for (i = 0; i < 26; i++){
		sha1hash12byte(str + i, prehash);

		hash[0] ^= prehash[0];
		hash[1] ^= prehash[1];
		hash[2] ^= prehash[2];
		hash[3] ^= prehash[3];
		hash[4] ^= prehash[4];
	}
}

#endif	// __SSE2__		//////////////////////////////////////////////////


int main(int argc, char *argv[])
{
	uint8_t data[100] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/ABCDEFGHIJKLMNOP";
	uint32_t hash[5] = {0};
	int i;

	sha1coinhash(data, hash);

	for (i = 0; i < 5; i++){
		printf("%02x%02x%02x%02x", (hash[i] >> 24) & 0xff, (hash[i] >> 16) & 0xff, (hash[i] >> 8) & 0xff, (hash[i]) & 0xff);
	}

	printf("\n");

	return 0;
}

