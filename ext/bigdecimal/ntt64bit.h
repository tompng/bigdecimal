#define NTT_BITS 64
#include "ntt.h"
#undef NTT_BITS

// TODO: uint128_t check
#define NTT64_PRIMITIVE_ROOT 5
#define NTT64_PRIME_BASE1 246
#define NTT64_PRIME_BASE2 247
#define NTT64_PRIME_SHIFT 56
#define NTT64_PRIME1 (((uint64_t)NTT64_PRIME_BASE1 << NTT64_PRIME_SHIFT) | 1)
#define NTT64_PRIME2 (((uint64_t)NTT64_PRIME_BASE2 << NTT64_PRIME_SHIFT) | 1)

// MAX_NTT64_BITS <= NTT64_PRIME_SHIFT && (1<<MAX_NTT64_BITS) * 99999999999999 ** 2 < (1<<128)
#define MAX_NTT64_BITS 34

/*
 * NTT multiplication
 * Uses two NTTs with mod (246 << 56 | 1), (247  << 56 | 1)
 * BASE: 10000000000000 (14 figures)
 */
void
ntt_multiply_14fig(uint64_t a_size, uint64_t b_size, uint64_t *a, uint64_t *b, uint64_t *c) {
    if (a_size < b_size) {
      ntt_multiply_14fig(b_size, a_size, b, a, c);
      return;
    }

    int b_bits = 0;
    while (((uint64_t)1 << b_bits) < b_size) b_bits++;
    int ntt_size_bits = b_bits + 1;
    if (ntt_size_bits > MAX_NTT64_BITS) {
      rb_raise(rb_eArgError, "Multiply size too large");
    }

    // To calculate large_a * small_b faster, split into several batches.
    uint64_t ntt_size = (uint64_t)1 << ntt_size_bits;
    uint64_t batch_size = ntt_size - (uint64_t)b_size;
    uint64_t batch_count = (uint64_t)((a_size + batch_size - 1) / batch_size);

    uint64_t *ntt1 = ruby_xcalloc(sizeof(uint64_t), ntt_size);
    uint64_t *ntt2 = ruby_xcalloc(sizeof(uint64_t), ntt_size);
    uint64_t *tmp1 = ruby_xcalloc(sizeof(uint64_t), ntt_size);
    uint64_t *tmp2 = ruby_xcalloc(sizeof(uint64_t), ntt_size);
    uint64_t *tmp3 = ruby_xcalloc(sizeof(uint64_t), ntt_size);
    uint64_t *conv1 = ruby_xcalloc(sizeof(uint64_t), ntt_size);
    uint64_t *conv2 = ruby_xcalloc(sizeof(uint64_t), ntt_size);

    // Calculate NTT for b in three primes
    // Result is reused for each batch of a.
    memcpy(tmp1, b, b_size * sizeof(uint64_t));
    memset(tmp1 + b_size, 0, (ntt_size - b_size) * sizeof(uint64_t));
    ntt64(ntt_size_bits, tmp1, ntt1, tmp2, NTT64_PRIMITIVE_ROOT, NTT64_PRIME_BASE1, NTT64_PRIME_SHIFT, +1);
    ntt64(ntt_size_bits, tmp1, ntt2, tmp2, NTT64_PRIMITIVE_ROOT, NTT64_PRIME_BASE2, NTT64_PRIME_SHIFT, +1);

    memset(c, 0, (a_size + b_size) * sizeof(uint64_t));
    for (uint64_t idx = 0; idx < batch_count; idx++) {
        if (idx == batch_count - 1) {
            uint64_t len = (uint64_t)a_size - idx * batch_size;
            memcpy(tmp1, a + idx * batch_size, len * sizeof(uint64_t));
            memset(tmp1 + len, 0, (ntt_size - len) * sizeof(uint64_t));
        } else {
            memcpy(tmp1, a + idx * batch_size, batch_size * sizeof(uint64_t));
        }
        // Calculate convolution for this batch in three primes
        ntt64(ntt_size_bits, tmp1, tmp2, tmp3, NTT64_PRIMITIVE_ROOT, NTT64_PRIME_BASE1, NTT64_PRIME_SHIFT, +1);
        for (uint64_t i = 0; i < ntt_size; i++) tmp2[i] = ((uint128_t)tmp2[i] * ntt1[i]) % NTT64_PRIME1;
        ntt64(ntt_size_bits, tmp2, conv1, tmp3, NTT64_PRIMITIVE_ROOT, NTT64_PRIME_BASE1, NTT64_PRIME_SHIFT, -1);
        ntt64(ntt_size_bits, tmp1, tmp2, tmp3, NTT64_PRIMITIVE_ROOT, NTT64_PRIME_BASE2, NTT64_PRIME_SHIFT, +1);
        for (uint64_t i = 0; i < ntt_size; i++) tmp2[i] = ((uint128_t)tmp2[i] * ntt2[i]) % NTT64_PRIME2;
        ntt64(ntt_size_bits, tmp2, conv2, tmp3, NTT64_PRIMITIVE_ROOT, NTT64_PRIME_BASE2, NTT64_PRIME_SHIFT, -1);

        // Restore the original convolution value from three convolutions calculated in three primes
        for (uint64_t i = 0; i < ntt_size; i++) {
            uint64_t mod1 = conv1[i], mod2 = conv2[i];
            // Calculate conv that satisfies: conv % PRIME1 == mod1 && conv % PRIME2 == mod2
            // conv = (mod1 * PRIME2 * PRIME2.pow(PRIME1 - 2, PRIME1) + mod2 * PRIME1 * PRIME1.pow(PRIME2 - 2, PRIME2)) % (PRIME1 * PRIME2)
            //      = (mod1 * (PRIME1 - 246) * PRIME2 + mod2 * (PRIME2 * 246 + 1)) % (PRIME1 * PRIME2)
            uint128_t conv = ((uint128_t)mod2 * 246 + (uint128_t)mod1 * (NTT64_PRIME1 - 246)) % NTT64_PRIME1 * NTT64_PRIME2 + mod2;
            uint64_t dl = conv % 100000000000000;
            uint64_t du = conv / 100000000000000;
            if (du) c[idx * batch_size + i] += du;
            if (dl) c[idx * batch_size + i + 1] += dl;
        }
    }
    uint64_t carry = 0;
    for (int64_t i = a_size + b_size - 1; i >= 0; i--) {
        uint64_t v = c[i] + carry;
        c[i] = v % 100000000000000;
        carry = v / 100000000000000;
    }
    ruby_xfree(ntt1);
    ruby_xfree(ntt2);
    ruby_xfree(tmp1);
    ruby_xfree(tmp2);
    ruby_xfree(tmp3);
    ruby_xfree(conv1);
    ruby_xfree(conv2);
}

const uint64_t bases14[14] = {
    1, 10, 100, 1000, 10000, 100000, 1000000, 10000000,
    100000000, 1000000000, 10000000000, 100000000000,
    1000000000000, 10000000000000
};
void
copy9to15(uint64_t *dst, uint32_t *src, uint64_t size) {
    for (uint64_t i = 0; i < size; i++) {
        uint32_t v = src[i];
        uint64_t start = i * 9;
        uint64_t end = start + 9;
        uint64_t js = start / 14;
        uint64_t je = end / 14;
        if (js == je) {
            dst[js] += v * bases14[5 - start % 14];
        } else {
            uint64_t b = bases14[start % 14 - 5];
            dst[js] += v / b;
            dst[je] += (v % b) * ((uint64_t)100000000000000 / b);
        }
    }
}

void
copy15to9(uint32_t *dst, uint64_t *src, uint64_t size) {
    for (uint64_t i = 0; i < size; i++) {
        uint64_t end = (i + 1) * 14 - 1;
        uint64_t j = end / 9;
        uint128_t v = (uint128_t)src[i] * bases14[8 - end % 9];
        uint32_t v0 = v % 1000000000;
        uint64_t v1 = v / 1000000000 % 1000000000;
        uint64_t v2 = v / 1000000000000000000;
        if (v2) dst[j - 2] += v2;
        if (v1) dst[j - 1] += v1;
        if (v0) dst[j] += v0;
    }
}

void
ntt_multiply64(uint64_t a_size, uint64_t b_size, uint32_t *a, uint32_t *b, uint32_t *c) {
    uint64_t a_size2 = (a_size * 9 + 13) / 14;
    uint64_t b_size2 = (b_size * 9 + 13) / 14;
    uint64_t c_size2 = a_size2 + b_size2;
    uint64_t *a2 = ruby_xcalloc(sizeof(uint64_t), a_size2);
    uint64_t *b2 = ruby_xcalloc(sizeof(uint64_t), b_size2);
    uint64_t *c2 = ruby_xcalloc(sizeof(uint64_t), c_size2);
    memset(a2, 0, sizeof(uint64_t) * a_size2);
    memset(b2, 0, sizeof(uint64_t) * b_size2);
    memset(c2, 0, sizeof(uint64_t) * c_size2);
    copy9to15(a2, a, a_size);
    copy9to15(b2, b, b_size);
    ntt_multiply_14fig(a_size2, b_size2, a2, b2, c2);
    memset(c, 0, sizeof(uint32_t) * (a_size + b_size));
    copy15to9(c, c2, c_size2);
    ruby_xfree(a2);
    ruby_xfree(b2);
    ruby_xfree(c2);
}
