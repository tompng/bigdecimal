#define NTT_BITS 32
#include "ntt.h"
#undef NTT_BITS

#define NTT32_PRIMITIVE_ROOT 17
#define NTT32_PRIME_BASE1 24
#define NTT32_PRIME_BASE2 26
#define NTT32_PRIME_BASE3 29
#define NTT32_PRIME_SHIFT 27
#define NTT32_PRIME1 (((uint32_t)NTT32_PRIME_BASE1 << NTT32_PRIME_SHIFT) | 1)
#define NTT32_PRIME2 (((uint32_t)NTT32_PRIME_BASE2 << NTT32_PRIME_SHIFT) | 1)
#define NTT32_PRIME3 (((uint32_t)NTT32_PRIME_BASE3 << NTT32_PRIME_SHIFT) | 1)
#define MAX_NTT32_BITS 27

/* Calculate c that satisfies: c % PRIME1 == mod1 && c % PRIME2 == mod2 && c % PRIME3 == mod3
 * Using Chinese Remainder Theorem, we can find c as follows:
 * c = (mod1 * 35002755423056150739595925972 + mod2 * 14584479687667766215746868453 + mod3 * 37919651490985126265126719818) % (PRIME1 * PRIME2 * PRIME3)
 */
static inline void
mod_restore_prime_24_26_29_shift_27(uint32_t mod1, uint32_t mod2, uint32_t mod3, uint32_t *digits) {
    // Use mixed radix notation to eliminate modulo by PRIME1 * PRIME2 * PRIME3
    // [DIG0, DIG1, DIG2] = DIG0 + DIG1 * PRIME1 + DIG2 * PRIME1 * PRIME2
    // DIG0: 0...PRIME1, DIG1: 0...PRIME2, DIG2: 0...PRIME3
    // 35002755423056150739595925972 = [1, 3489660916, 3113851359]
    // 14584479687667766215746868453 = [0, 13, 1297437912]
    // 37919651490985126265126719818 = [0, 0, 3373338954]
    uint64_t c0 = mod1;
    uint64_t c1 = (uint64_t)mod2 * 13 + (uint64_t)mod1 * 3489660916;
    uint64_t c2 = (uint64_t)mod3 * 3373338954 % NTT32_PRIME3 + (uint64_t)mod2 * 1297437912 % NTT32_PRIME3 + (uint64_t)mod1 * 3113851359 % NTT32_PRIME3;
    c2 += c1 / NTT32_PRIME2;
    c1 %= NTT32_PRIME2;
    c2 %= NTT32_PRIME3;
    // Base conversion
    uint32_t base = 1000000000;
    c1 += c2 % base * NTT32_PRIME2;
    c0 += c1 % base * NTT32_PRIME1;
    c1 /= base;
    digits[0] = c0 % base;
    c0 /= base;
    c1 += c2 / base % base * NTT32_PRIME2;
    c0 += c1 % base * NTT32_PRIME1;
    c1 /= base;
    digits[1] = c0 % base;
    c0 = c0 / base + c1 % base * NTT32_PRIME1;
    digits[2] = c0 % base;
    digits[3] = (c0 / base + c1 / base % base * NTT32_PRIME1) % base;
}

/*
 * NTT multiplication
 * Uses three NTTs with mod (24 << 27 | 1), (26 << 27 | 1), and (29 << 27 | 1)
 */
static void
ntt_multiply32(uint32_t a_size, uint32_t b_size, uint32_t *a, uint32_t *b, uint32_t *c) {
    if (a_size < b_size) {
      ntt_multiply32(b_size, a_size, b, a, c);
      return;
    }

    int b_bits = 0;
    while (((uint32_t)1 << b_bits) < b_size) b_bits++;
    int ntt_size_bits = b_bits + 1;
    if (ntt_size_bits > MAX_NTT32_BITS) {
      rb_raise(rb_eArgError, "Multiply size too large");
    }

    // To calculate large_a * small_b faster, split into several batches.
    uint32_t ntt_size = (uint32_t)1 << ntt_size_bits;
    uint32_t batch_size = ntt_size - (uint32_t)b_size;
    uint32_t batch_count = (uint32_t)((a_size + batch_size - 1) / batch_size);

    uint32_t *ntt1 = ruby_xcalloc(sizeof(uint32_t), ntt_size);
    uint32_t *ntt2 = ruby_xcalloc(sizeof(uint32_t), ntt_size);
    uint32_t *ntt3 = ruby_xcalloc(sizeof(uint32_t), ntt_size);
    uint32_t *tmp1 = ruby_xcalloc(sizeof(uint32_t), ntt_size);
    uint32_t *tmp2 = ruby_xcalloc(sizeof(uint32_t), ntt_size);
    uint32_t *tmp3 = ruby_xcalloc(sizeof(uint32_t), ntt_size);
    uint32_t *conv1 = ruby_xcalloc(sizeof(uint32_t), ntt_size);
    uint32_t *conv2 = ruby_xcalloc(sizeof(uint32_t), ntt_size);
    uint32_t *conv3 = ruby_xcalloc(sizeof(uint32_t), ntt_size);

    // Calculate NTT for b in three primes
    // Result is reused for each batch of a.
    memcpy(tmp1, b, b_size * sizeof(uint32_t));
    memset(tmp1 + b_size, 0, (ntt_size - b_size) * sizeof(uint32_t));
    ntt32(ntt_size_bits, tmp1, ntt1, tmp2, NTT32_PRIMITIVE_ROOT, NTT32_PRIME_BASE1, NTT32_PRIME_SHIFT, +1);
    ntt32(ntt_size_bits, tmp1, ntt2, tmp2, NTT32_PRIMITIVE_ROOT, NTT32_PRIME_BASE2, NTT32_PRIME_SHIFT, +1);
    ntt32(ntt_size_bits, tmp1, ntt3, tmp2, NTT32_PRIMITIVE_ROOT, NTT32_PRIME_BASE3, NTT32_PRIME_SHIFT, +1);

    memset(c, 0, (a_size + b_size) * sizeof(uint32_t));
    for (uint32_t idx = 0; idx < batch_count; idx++) {
        if (idx == batch_count - 1) {
            uint32_t len = (uint32_t)a_size - idx * batch_size;
            memcpy(tmp1, a + idx * batch_size, len * sizeof(uint32_t));
            memset(tmp1 + len, 0, (ntt_size - len) * sizeof(uint32_t));
        } else {
            memcpy(tmp1, a + idx * batch_size, batch_size * sizeof(uint32_t));
        }
        // Calculate convolution for this batch in three primes
        ntt32(ntt_size_bits, tmp1, tmp2, tmp3, NTT32_PRIMITIVE_ROOT, NTT32_PRIME_BASE1, NTT32_PRIME_SHIFT, +1);
        for (uint32_t i = 0; i < ntt_size; i++) tmp2[i] = ((uint64_t)tmp2[i] * ntt1[i]) % NTT32_PRIME1;
        ntt32(ntt_size_bits, tmp2, conv1, tmp3, NTT32_PRIMITIVE_ROOT, NTT32_PRIME_BASE1, NTT32_PRIME_SHIFT, -1);
        ntt32(ntt_size_bits, tmp1, tmp2, tmp3, NTT32_PRIMITIVE_ROOT, NTT32_PRIME_BASE2, NTT32_PRIME_SHIFT, +1);
        for (uint32_t i = 0; i < ntt_size; i++) tmp2[i] = ((uint64_t)tmp2[i] * ntt2[i]) % NTT32_PRIME2;
        ntt32(ntt_size_bits, tmp2, conv2, tmp3, NTT32_PRIMITIVE_ROOT, NTT32_PRIME_BASE2, NTT32_PRIME_SHIFT, -1);
        ntt32(ntt_size_bits, tmp1, tmp2, tmp3, NTT32_PRIMITIVE_ROOT, NTT32_PRIME_BASE3, NTT32_PRIME_SHIFT, +1);
        for (uint32_t i = 0; i < ntt_size; i++) tmp2[i] = ((uint64_t)tmp2[i] * ntt3[i]) % NTT32_PRIME3;
        ntt32(ntt_size_bits, tmp2, conv3, tmp3, NTT32_PRIMITIVE_ROOT, NTT32_PRIME_BASE3, NTT32_PRIME_SHIFT, -1);

        // Restore the original convolution value from three convolutions calculated in three primes
        for (uint32_t i = 0; i < ntt_size; i++) {
            uint32_t dig[4];
            mod_restore_prime_24_26_29_shift_27(conv1[i], conv2[i], conv3[i], dig);
            for (int j = 0; j < 4; j++) {
                // Maximum overlap(4) * maximum value(999999999) does not overflow 32-bit integer
                // Index check: if dig[j] is non-zero, assign index is valid
                if (dig[j]) c[idx * batch_size + i + 1 - j] += dig[j];
            }
        }
    }
    uint32_t carry = 0;
    for (int32_t i = a_size + b_size - 1; i >= 0; i--) {
        uint32_t v = c[i] + carry;
        c[i] = v % 1000000000;
        carry = v / 1000000000;
    }
    ruby_xfree(ntt1);
    ruby_xfree(ntt2);
    ruby_xfree(ntt3);
    ruby_xfree(tmp1);
    ruby_xfree(tmp2);
    ruby_xfree(tmp3);
    ruby_xfree(conv1);
    ruby_xfree(conv2);
    ruby_xfree(conv3);
}

static void
ntt_multiply(size_t a_size, size_t b_size, uint32_t *a, uint32_t *b, uint32_t *c) {
    size_t size = Min(a_size, b_size);
    if (size <= (size_t)1 << (MAX_NTT32_BITS - 1)) {
        ntt_multiply32((uint32_t)a_size, (uint32_t)b_size, a, b, c);
    } else {
        ntt_multiply64((uint64_t)a_size, (uint64_t)b_size, a, b, c);
    }
}
