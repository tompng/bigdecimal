#define MOD_POW NTT_FUNCNAME(mod_pow)
#define NTT_RECURSIVE NTT_FUNCNAME(ntt_recursive)
#define NTT NTT_FUNCNAME(ntt)

// Calculates base**ex % mod
static UINT
MOD_POW(UINT base, UINT ex, UINT mod) {
    UINT res = 1;
    UINT bit = 1;
    while (true) {
        if (ex & bit) {
            ex ^= bit;
            res = MULT_MOD(res, base, mod);
        }
        if (!ex) break;
        base = MULT_MOD(base, base, mod);
        bit <<= 1;
    }
    return res;
}

// Recursively performs butterfly operations of NTT
static void
NTT_RECURSIVE(int size_bits, UINT *input, UINT *output, UINT *tmp, int depth, UINT r, UINT prime) {
    if (depth > 0) {
        NTT_RECURSIVE(size_bits, input, tmp, output, depth - 1, MULT_MOD(r, r, prime), prime);
    } else {
        tmp = input;
    }
    UINT size_half = (UINT)1 << (size_bits - 1);
    UINT stride = (UINT)1 << (size_bits - depth - 1);
    UINT n = size_half / stride;
    UINT rn = 1, rm = prime - 1;
    UINT idx = 0;
    for (UINT i = 0; i < n; i++) {
        UINT j = i * 2 * stride;
        for (UINT k = 0; k < stride; k++, j++, idx++) {
            UINT a = tmp[j], b = tmp[j + stride];
            output[idx] = ADD_MULT_MOD(a, rn, b, prime);
            output[idx + size_half] = ADD_MULT_MOD(a, rm, b, prime);
        }
        rn = MULT_MOD(rn, r, prime);
        rm = MULT_MOD(rm, r, prime);
    }
}

/* Perform NTT on input array.
 * base, shift: Represent the prime number as (base << shift | 1)
 * r_base: Primitive root of unity modulo prime
 * size_bits: log2 of the size of the input array. Should be less or equal to shift
 * input: input array of size 1 << size_bits
 */
static void
NTT(int size_bits, UINT *input, UINT *output, UINT *tmp, int r_base, int base, int shift, int dir) {
    UINT size = (UINT)1 << size_bits;
    UINT prime = ((UINT)base << shift) | 1;

    // rmax**(1 << shift) % prime == 1
    // r**size % prime == 1
    UINT rmax = MOD_POW(r_base, base, prime);
    UINT r = MOD_POW(rmax, (UINT)1 << (shift - size_bits), prime);

    if (dir < 0) r = MOD_POW(r, prime - 2, prime);
    NTT_RECURSIVE(size_bits, input, output, tmp, size_bits - 1, r, prime);
    if (dir < 0) {
        UINT n_inv = MOD_POW((UINT)size, prime - 2, prime);
        for (UINT i = 0; i < size; i++) {
            output[i] = MULT_MOD(output[i], n_inv, prime);
        }
    }
}

#undef MOD_POW
#undef NTT_RECURSIVE
#undef NTT
