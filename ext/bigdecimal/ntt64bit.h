#define NTT64_PRIMITIVE_ROOT 3
#define NTT64_PRIME_BASE 5
#define NTT64_PRIME_SHIFT 55
#define NTT64_PRIME ((uint64_t)NTT64_PRIME_BASE << NTT64_PRIME_SHIFT | 1)
#define NTT64_DECDIG_BASE 1000
// Constraints: (NTT64_DECDIG_BASE-1)**2*(1<<MAX_NTT64_BITS) < 1<<NTT64_PRIME_SHIFT
#define MAX_NTT64_BITS 35

static inline uint64_t
add_mod_ntt64_prime(uint64_t a, uint64_t b) {
    return NTT64_PRIME - a <= b ? a + b - NTT64_PRIME : a + b;
}

static inline uint64_t
mult_mod_ntt64_prime(uint64_t a, uint64_t b) {
    uint32_t a1 = a >> 32, a0 = a & 0xffffffff;
    uint32_t b1 = b >> 32, b0 = b & 0xffffffff;
    uint64_t ab00 = (uint64_t)a0 * b0;
    uint64_t ab01 = (uint64_t)a0 * b1;
    uint64_t ab10 = (uint64_t)a1 * b0;
    uint64_t ab11 = (uint64_t)a1 * b1;
    uint64_t tmp = (ab00 >> 32) + (ab01 & 0xffffffff)+ (ab10 & 0xffffffff);
    uint64_t ab0 = (ab00 & 0xffffffff) | (tmp & 0xffffffff) << 32;
    uint64_t ab1 = (tmp >> 32) + (ab01 >> 32) + (ab10 >> 32) + ab11;
    // Maximum value of ab1 is (NTT64_PRIME_BASE**2) << (2 * NTT64_PRIME_SHIFT - 64)
    // ab1 << (64 - NTT64_PRIME_SHIFT) does not overflow uint64_t when NTT64_PRIME_BASE==5 and NTT64_PRIME_SHIFT==55
    uint64_t div = (ab1 << (64 - NTT64_PRIME_SHIFT) | ab0 >> NTT64_PRIME_SHIFT) / NTT64_PRIME_BASE;
    uint64_t mod = div ? ab0 - (div - 1) * NTT64_PRIME : ab0;
    return mod >= NTT64_PRIME ? mod - NTT64_PRIME : mod;
}

#define NTT_FUNCNAME(name) name##64
#define UINT uint64_t
#define ADD_MULT_MOD(a, b, c, mod) add_mod_ntt64_prime(a, mult_mod_ntt64_prime(b, c))
#define MULT_MOD(a, b, mod) mult_mod_ntt64_prime(a, b)
#include "ntt.h"
#undef UINT
#undef NTT_FUNCNAME
#undef ADD_MULT_MOD
#undef MULT_MOD

void
ntt_multiply64_3fig(uint64_t a_size, uint64_t b_size, uint64_t *a, uint64_t *b, uint64_t *c) {
    if (a_size < b_size) {
      ntt_multiply64_3fig(b_size, a_size, b, a, c);
      return;
    }
    size_t c_size = a_size + b_size;

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

    uint64_t *ntt = ruby_xcalloc(sizeof(uint64_t), ntt_size);
    uint64_t *tmp1 = ruby_xcalloc(sizeof(uint64_t), ntt_size);
    uint64_t *tmp2 = ruby_xcalloc(sizeof(uint64_t), ntt_size);
    uint64_t *conv = ruby_xcalloc(sizeof(uint64_t), ntt_size);

    // Calculate NTT for b in three primes
    // Result is reused for each batch of a.
    memcpy(tmp1, b, b_size * sizeof(uint64_t));
    memset(tmp1 + b_size, 0, (ntt_size - b_size) * sizeof(uint64_t));
    ntt64(ntt_size_bits, tmp1, ntt, tmp2, NTT64_PRIMITIVE_ROOT, NTT64_PRIME_BASE, NTT64_PRIME_SHIFT, +1);

    memset(c, 0, c_size * sizeof(uint64_t));
    for (uint64_t idx = 0; idx < batch_count; idx++) {
        uint64_t len = idx == batch_count - 1 ? a_size - idx * batch_size : batch_size;
        memcpy(tmp1, a + idx * batch_size, len * sizeof(uint64_t));
        memset(tmp1 + len, 0, (ntt_size - len) * sizeof(uint64_t));
        // Calculate convolution for this batch in three primes
        ntt64(ntt_size_bits, tmp1, tmp2, conv, NTT64_PRIMITIVE_ROOT, NTT64_PRIME_BASE, NTT64_PRIME_SHIFT, +1);
        for (uint64_t i = 0; i < ntt_size; i++) tmp1[i] = mult_mod_ntt64_prime(tmp2[i], ntt[i]);
        ntt64(ntt_size_bits, tmp1, conv, tmp2, NTT64_PRIMITIVE_ROOT, NTT64_PRIME_BASE, NTT64_PRIME_SHIFT, -1);
        for (uint64_t i = 0, j = idx * batch_size + 1; i < ntt_size; i++, j++) {
            if (j < c_size) c[j] += conv[i];
        }
    }
    uint64_t carry = 0;
    for (int64_t i = c_size - 1; i >= 0; i--) {
        uint64_t v = c[i] + carry;
        c[i] = v % NTT64_DECDIG_BASE;
        carry = v / NTT64_DECDIG_BASE;
    }
    ruby_xfree(ntt);
    ruby_xfree(tmp1);
    ruby_xfree(tmp2);
    ruby_xfree(conv);
}

void
ntt_multiply64(uint64_t a_size, uint64_t b_size, uint32_t *a, uint32_t *b, uint32_t *c) {
    uint64_t a_size2 = a_size * 3;
    uint64_t b_size2 = b_size * 3;
    uint64_t c_size = a_size + b_size;
    uint64_t c_size2 = c_size * 3;
    uint64_t *a2 = ruby_xcalloc(sizeof(uint64_t), a_size2);
    uint64_t *b2 = ruby_xcalloc(sizeof(uint64_t), b_size2);
    uint64_t *c2 = ruby_xcalloc(sizeof(uint64_t), c_size2);
    for (uint64_t i = 0; i < a_size; i++) {
        a2[i * 3] = a[i] / 1000000;
        a2[i * 3 + 1] = (a[i] / 1000) % 1000;
        a2[i * 3 + 2] = a[i] % 1000;
    }
    for (uint64_t i = 0; i < b_size; i++) {
        b2[i * 3] = b[i] / 1000000;
        b2[i * 3 + 1] = (b[i] / 1000) % 1000;
        b2[i * 3 + 2] = b[i] % 1000;
    }
    memset(c2, 0, sizeof(uint64_t) * c_size2);
    ntt_multiply64_3fig(a_size2, b_size2, a2, b2, c2);
    for (uint64_t i = 0; i < c_size; i++) {
        c[i] = (uint32_t)(c2[i * 3] * 1000000 + c2[i * 3 + 1] * 1000 + c2[i * 3 + 2]);
    }
    ruby_xfree(a2);
    ruby_xfree(b2);
    ruby_xfree(c2);
}
