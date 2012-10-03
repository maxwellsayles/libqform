/**
 * 64bit negative discriminant qforms
 *
 * compose/square can handle 59 bit discriminants
 * cube can handle 59 bit discriminants
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gcdext_binary_l2r.h"
#include "math32.h"
#include "math64.h"
#include "math_mpz.h"
#include "mpz_qform.h"
#include "primes.h"
#include "s64_qform.h"
#include "sqrtmodp_list.h"

const group_cost_t s64_qform_costs = {
    1.218,
    1.000,
    1.720
};

#define gcdext_s32(s, t, u, v) gcdext_binary_l2r_s32(s, t, u, v)
#define gcdext_left_s32(s, u, v) gcdext_binary_l2r_s32(s, 0, u, v)
#define gcdext_partial_s32(R1, R0, C1, C0, bound) gcdext_partial_binary_l2r_s32(R1, R0, C1, C0, bound)
#define gcdext_s64(s, t, u, v) gcdext_binary_l2r_s64(s, t, u, v)
#define gcdext_left_s64(s, u, v) gcdext_binary_l2r_s64(s, 0, u, v)
#define gcdext_partial_s64(R1, R0, C1, C0, bound) gcdext_partial_binary_l2r_s64(R1, R0, C1, C0, bound)

/**
 * initialize the group
 */
void s64_qform_group_init(s64_qform_group_t* group) {
    group->desc.group.elem_init = (group_elem_init_f*)&s64_qform_init;
    group->desc.group.elem_clear = (group_elem_clear_f*)&s64_qform_clear;
    group->desc.group.elem_size = sizeof(s64_qform_t);

    group->desc.group.hash32 = (group_hash32_f*)&s64_qform_hash32;
    group->desc.group.set_id = (group_set_id_f*)&s64_qform_set_id;
    group->desc.group.is_id = (group_is_id_f*)&s64_qform_is_id;
    group->desc.group.set = (group_set_f*)&s64_qform_set;
    group->desc.group.equal = (group_equal_f*)&s64_qform_equal;
    group->desc.group.inverse = (group_inverse_f*)&s64_qform_inverse;
    group->desc.group.compose = (group_compose_f*)&s64_qform_compose;
    group->desc.group.square = (group_square_f*)&s64_qform_square;
    group->desc.group.cube = (group_cube_f*)&s64_qform_cube;
    group->desc.group.print = (group_print_f*)&s64_qform_print;

    group->desc.discriminant_max_bits = s64_qform_group_max_bits;
    group->desc.clear = (qform_group_clear_f*)&s64_qform_group_clear;
    group->desc.set_discriminant = (qform_group_set_discriminant_f*)&s64_qform_group_set_discriminant;
    group->desc.reduce = (qform_reduce_f*)&s64_qform_reduce;
    group->desc.is_primeform = (qform_is_primeform_f*)&s64_qform_is_primeform;
    group->desc.is_ambiguous = (qform_is_ambiguous_f*)&s64_qform_is_ambiguous;
    group->desc.split_ambiguous = (qform_split_ambiguous_f*)&s64_qform_split_ambiguous;
}



static inline int32_t avg_s32(const int32_t a, const int32_t b) {
#if defined(__x86_64) || (defined(__i386) && !defined(__APPLE__))
    int32_t res;
    asm(
        "movl %%eax, %%edx\n\t"
        "movl %%ebx, %%ecx\n\t"
        "sarl $31, %%edx\n\t"
        "sarl $31, %%ecx\n\t"
        "addl %%ebx, %%eax\n\t"
        "adcl %%ecx, %%edx\n\t"
        "shrdl $1, %%edx, %%eax\n\t"
        : "=&a"(res)
        : "0"(a), "b"(b)
        : "ecx", "edx", "cc"
    );
    return res;
#else
    return ((int64_t)a + (int64_t)b) >> 1;
#endif
}



// res=(f1*f2+f3*f4)/d
static inline int64_t muladdmuldiv_s64(const int64_t f1, const int64_t f2, const int64_t f3, const int64_t f4, const int64_t d) {
    int64_t res;
#if defined(__x86_64)
    asm(
        "movq %1, %%rax\n"
        "imulq %2\n"
        "movq %%rax, %%r10\n"
        "movq %3, %%rax\n"
        "movq %%rdx, %%r11\n"
        "imulq %4\n"
        "addq %%r10, %%rax\n"
        "adcq %%r11, %%rdx\n"
        "idivq %5\n"
        : "=&a"(res)
        : "rm"(f1), "rm"(f2), "rm"(f3), "rm"(f4), "rm"(d)
        : "rdx", "r10", "r11", "cc"
    );
#else
    s128_t t1;
    s128_t t2;
    s128_t q;
    mul_s128_s64_s64(&t1, f1, f2);
    mul_s128_s64_s64(&t2, f3, f4);
    add_s128_s128(&t1, &t2);
    div_s128_s128_s64(&q, &t1, d);
    res = get_s64_from_s128(&q);
#endif
    return res;
}


/**
 * r = a >> (size(a)/2) and a > 0
 * r ~  SquareRoot(a)
 */
static inline uint32_t half_rshift_u32(uint32_t a) {
    uint32_t r = a;
#if defined(__i386) || defined(__x86_64)
    asm(
        "xorl %%ecx, %%ecx\n\t"
        "bsrl %0, %%ecx\n\t"
        "addl $1, %%ecx\n\t"
        "shrl $1, %%ecx\n\t"
        "shrl %%cl, %0"
        : "=rm"(r)
        : "0"(r)
        : "cc", "ecx"
    );
#else
    r >>= ((msb_u32(r)+1)>>1;
#endif
    // should not return 0 since this causes weirdness with partial gcd
    if (r == 0) return 1;
    return r;
}


// c = (b^2-D)/(4a)
static inline int64_t s64_qform_c(const s64_qform_group_t* group, int32_t a, int32_t b) {
    return (((int64_t)b*(int64_t)b - group->D) >> 2) / a;
}



void s64_qform_set_id(s64_qform_group_t* group, s64_qform_t* form) {
    form->a = 1;
    if ((group->D & 3) == 0) {
        // D is divisible by four
        form->b = 0;
    }
    else {
        form->b = 1;
    }

    // c = (b*b-D) / (4a)
    form->c = s64_qform_c(group, form->a, form->b);
}


/**
 * Check for a primeform where a=p.
 * Tests b = +/- sqrt(D) mod p.
 * Note: -0 mod p is p.  This way we try ambiguous forms.
 * @return 0 if there is no prime form and 1 otherwise.
 */
int s64_qform_is_primeform(s64_qform_group_t* group, s64_qform_t* form, const int p) {
    int Dmodp;
    const short* sqrtp;

    if (p > sqrtmodp_maxp)
        return 0; // p is too large for the table

    sqrtp = sqrtmodp[p];
    if (sqrtp == 0)
        return 0; // p was not prime

    // a = p
    form->a = p;

    // b = sqrt(D) (mod p)
    Dmodp = group->D % p;
    if (Dmodp < 0) {
        Dmodp += p;
    }
    if (Dmodp == 0) {
        // prime divides the discriminant
        return 0;
    }

    form->b = sqrtp[Dmodp];
    if (form->b == -1) {
        return 0; // D does not have a sqrt (mod p)
    }    

    // We know that p | b^2-D for +/- b.
    // Compute c = (b^2 - d) / (4a) if possible.
    if (p == 2) {
        // special case if p == 2
        form->c = (int64_t)form->b * (int64_t)form->b;
        form->c -= group->D;
        if ((form->c & 7) == 0) {
            // 4a | b^2-D, a=2, so we divide by 8
            form->c >>= 3;
            return 1;
        }
        if (form->b == 1)
            return 0;

        // try -b = 2
        form->b = 2;
        form->c = 4 - group->D;
        if ((form->c & 7) != 0)
            return 0;

        // 4a | b^2-D, a=2, so we divide by 8
        form->c >>= 3;
        return 1;
    }
    // p != 2

    form->c = (int64_t)form->b * (int64_t)form->b;
    form->c -= group->D;

    if ((form->c & 3) == 0) {
        // 4a | b^2-D

        // divide by 4
        form->c >>= 2;

        // divide by a
        form->c /= p;
        return 1;
    }

    // try -b
    form->b = p-form->b;

    // make sure that 4 | b^2-D
    form->c = (int64_t)form->b * (int64_t)form->b;
    form->c -= group->D;
    if ((form->c & 3) != 0)
        return 0;

    // divide by 4
    form->c >>= 2;

    // divide by a
    form->c /= p;
    
    return 1;
}


/**
 * If N is split by the form (a,b,C), then d is a factor and 1 is returned
 * otherwise 0 is returned.
 * Assumes that form is an ambiguous form.
 */
int s64_qform_split_ambiguous(s64_qform_group_t* group, mpz_t out_d, const mpz_t in_N, const s64_qform_t* form) {
    int64_t m;
    int64_t N = mpz_get_s64(in_N);
    int64_t d;

    // ambiguous forms are of three types
    // (a, 0, c), (a, a, c), and (a, b, a)
    if (form->b == 0 || form->a == form->b) {
        // first two forms (a, 0, c) and (a, a, c)
        // (a, 0, c) => -4ac = D => a | D
        // (a, a, c) => a^2-4ac = a(a-4c) = D => a | d
        m = form->a;
    } else {
        // final form (a, b, a)
        // (a, b, a) => b^2 - 4a^2 = (b+2a)(b-2a) = D
        // pick m = min |b +/- 2a|
        if (form->b < 0) {
            m = form->b + 2*form->a;
        } else {
            m = form->b - 2*form->a;
        }
    }

    if (m < 0) m = -m;
    d = gcd_binary_l2r_u64(m, N);
    mpz_set_s64(out_d, d);

    // return true if this is a non-trivial root
    return (d > 1 && d < N);
}


/**
 * Saves the discriminant and computes the fourth root.
 * Assumes that D is negative.
 */
void s64_qform_group_set_discriminant(s64_qform_group_t* group, const mpz_t in_D) {
    s64_qform_group_set_discriminant_s64(group, mpz_get_s64(in_D));
}

void s64_qform_group_set_discriminant_s64(s64_qform_group_t* group, const int64_t D) {
    group->D = D;
    group->S = sqrt_u64(abs_s64(D));
    group->L = sqrt_u64(group->S);
}


// a form is reduced when -a < b <= a < c or 0 <= b <= a = c
void s64_qform_reduce(s64_qform_group_t* group, s64_qform_t* form) {
    int32_t q;
    int32_t r;
    int looping = 1;

    while (looping) {
        looping = 0;
        if (form->a > form->c) {
            // swap a and c and invert b
            r = form->a;
            form->a = form->c;
            form->c = r;
            form->b = -form->b;
        }
        if ((form->b > form->a) || (form->b <= -form->a)) {

            // find r such that -a < r <= a
            // and r = b (mod 2a).
            // q = b/2a = (b/a)/2
            // r = b%2a = (b%a) +- a*((b/a)%2)
            // where '/' is divide with floor
            divrem_s32(&q, &r, form->b, form->a);
            if (q&1) {
                if (r <= 0) {
                    // r is negative
                    if (form->a > 0) {
                        r += form->a;
                        -- q;
                    }
                    else {
                        r -= form->a;
                        ++ q;
                    }
                }
                else {
                    // r is positive
                    if (form->a > 0) {
                        r -= form->a;
                        ++ q;
                    }
                    else {
                        r += form->a;
                        -- q;
                    }
                }
            }
            q >>= 1;

            form->c -= ((int64_t)q * ((int64_t)form->b + (int64_t)r)) >> 1;
            
            form->b = r;
            
            looping = 1;
        }
    }

    // account for special case
    if ((form->a == form->c) && (form->b < 0)) {
        form->b = -form->b;
    }
}


/**
 * nucomp algorithm from "Solving the Pell Equation"
 */
void s64_qform_compose(s64_qform_group_t* group, s64_qform_t* C, const s64_qform_t* A, const s64_qform_t* B) {
    int32_t a1, a2, b1, b2;
    int64_t c2;
    int32_t g, s, x, y, z;
    int32_t m12, p12, u;
    int32_t C1, C0, r1, r0;
    int32_t m1, m2;
    int64_t bound;
    int64_t tmp64;

    // if A is the identity, return B
    if (s64_qform_is_id(group, A)) {
        s64_qform_set(group, C, B);
        s64_qform_reduce(group, C);
        return;
    }

    // if B is the identity, return A
    if (s64_qform_is_id(group, B)) {
        s64_qform_set(group, C, A);
        s64_qform_reduce(group, C);
        return;
    }

    // Local copies
    if (A->a > B->a) {
        a1 = A->a;
        b1 = A->b;
        a2 = B->a;
        b2 = B->b;
        c2 = B->c;
    } else {
        a1 = B->a;
        b1 = B->b;
        a2 = A->a;
        b2 = A->b;
        c2 = A->c;
    }

    // Compute d1=gcd(a1, a2) and u1 such that  $a2 | (d1 - a2 * u1)$
    g = gcdext_left_s32(&x, a2, a1); 

    // Compute gcd((b1+b2)/2, g) = s = y * (b1+b2)/2 + z * g
    p12 = avg_s32(b1, b2);
    m12 = p12 - b2;

    // Compute u = x*z*(b1-b2)/2 - y*c mod a1
    u = mulmod_s32(x, m12, a1);
    s = 1;
    if (g != 1) {
        s = gcdext_s32(&y, &z, p12, g);
        if (s != 1) {
            a1 /= s;
            a2 /= s;
        }

        u = mulmod_s32(u, z, a1);
        r0 = c2 % a1;
        r1 = mulmod_s32(y, r0, a1);
        u = submod_s32(u, r1, a1);
    }
    if (u < 0) {
        u += a1;
    }

    // compute bounds
    bound = (half_rshift_u32(a1/a2)) * group->L;

    if (a1 <= bound) {
        // normal composition
        C->a = a1*a2;

        // C->b = 2*a2*u + b2 (mod 2C->a)
        tmp64 = ((int64_t)a2*(int64_t)u << 1) + b2;
        C->b = mod_s32_s64_u32(tmp64, C->a << 1);
    } 
    else {
        // nucomp steps
        r1 = a1;
        r0 = u;

        // partial xgcd
        gcdext_partial_s32((uint32_t*)&r1, (uint32_t*)&r0, &C1, &C0, bound);

        // m1 = (a2*r0 + m12*C0) / a1
        m1 = muladdmul_s64_4s32(a2, r0, m12, C0) / a1;

        // m2 = (p12*r0 - s*C0*c2) / a1
        m2 = ((int64_t)p12 * (int64_t)r0 - (int64_t)s * (int64_t)C0 * c2) / a1;
        
        // a_{i+1} = r0*m1 - C0*m2
        C->a = r0*m1 - C0*m2;

        // b_{i+1} = (a2*r0 - a*|C1|)/C0 (mod 2a)
        if (C1 >= 0) {
            C1 = -C1;
        }
        tmp64 = muladdmul_s64_4s32(a2, r0, C->a, C1) / C0;
        tmp64 = (tmp64<<1)-b2;
        C->b = mod_s32_s64_u32(tmp64, C->a << 1);
    }

    C->c = s64_qform_c(group, C->a, C->b);
    s64_qform_reduce(group, C);
}


/**
 * nudupl, adapted from compose above
 */
void s64_qform_square(s64_qform_group_t* group, s64_qform_t* C, const s64_qform_t* A) {
    int32_t a1, b1;
    int64_t c1;
    int32_t s, y;
    int32_t u;
    int32_t C1, C0, r1, r0, m2;
    int64_t tmp64;

    // if A is the identity, return the identity
    if (s64_qform_is_id(group, A)) {
        s64_qform_set_id(group, C);
        return;
    }

    // Local copies
    a1 = A->a;
    b1 = A->b;
    c1 = A->c;

    // Compute d1=gcd(a1, a2) and u1 such that  $a2 | (d1 - a2 * u1)$

    // Compute gcd((b1+b2)/2, g) = s = y * (b1+b2)/2 + z * g

    // Compute u = x*z*(b1-b2)/2 - y*c mod a1
    s = gcdext_left_s32(&y, b1, a1);
    a1 /= s;

    r0 = c1 % a1;
    r1 = mulmod_s32(y, r0, a1);
    u = -r1;
    if (u < 0) {
        u += a1;
    }

    if (a1 <= group->L) {
        // normal composition
        C->a = a1*a1;

        // C->b = 2*a2*u + b2 (mod 2a)
        tmp64 = ((int64_t)a1*(int64_t)u << 1) + b1;
        C->b = mod_s32_s64_u32(tmp64, C->a << 1);
    }
    else {
        // nucomp steps

        // partial xgcd
        r1 = a1;
        r0 = u;
        gcdext_partial_s32((uint32_t*)&r1, (uint32_t*)&r0, &C1, &C0, group->L);

        // m2 = (b1 * r0 - s*C0*c1) / a1
        m2 = ((int64_t)b1 * (int64_t)r0 - (int64_t)s * (int64_t)C0 * c1) / a1;

        // a_{i+1} = r0^2 - C0*m2
        C->a = r0*r0 - C0*m2;

        // b_{i+1} = (a1 * r0 - a * |C1|)/C0  (mod 2a)
        if (C1 >= 0) {
            C1 = -C1;
        }
        tmp64 = muladdmul_s64_4s32(a1, r0, C->a, C1) / C0;
        tmp64 = (tmp64<<1) - b1;
        C->b = mod_s32_s64_u32(tmp64, C->a << 1);
    }

    C->c = s64_qform_c(group, C->a, C->b);
    s64_qform_reduce(group, C);
}



/**
 * Computes a reduced ideal equivalent to the cube of an ideal
 * using an adaptation of Shanks' NUCOMP algorithm.
 */
void s64_qform_cube(s64_qform_group_t* group, s64_qform_t* R, const s64_qform_t* A) {

    int32_t a1;
    int32_t b1;
    int64_t c1;

    int32_t SP;
    int32_t v1;
    int32_t N;
    int32_t T;
    int32_t S, u2, v2;
    int32_t R1, C1, C2;
    int32_t M1, M2;

    int32_t L_32;
    int32_t t1_32, t2_32;

    int64_t temp, temp2;

    int64_t B;
    int64_t L;
    int64_t K;
    int64_t u2_64, v2_64;
    int64_t R1_64, R2_64;
    int64_t C1_64, C2_64;
    int64_t T_64;


    // if A is the identity, return the identity
    if (s64_qform_is_id(group, A)) {
        s64_qform_set_id(group, R);
        return;
    }

    // cubing an ambiguous form gets us the form back
    if (s64_qform_is_ambiguous(group, A)) {
        s64_qform_set(group, R, A);
        s64_qform_reduce(group, R);
        return;
    }


    a1 = A->a;
    b1 = A->b;
    c1 = A->c;

    // solve SP = v1 b + u1 a (only need v1)
    SP = gcdext_left_s32(&v1, b1, a1);

    if (SP == 1) {

        // N = a
        N = a1;

        // L = a^2
        L = (int64_t)a1 * (int64_t)a1;

        // K = c v1 (v1(b - a c v1) - 2) (mod L)
        if (s64_is_s32(L)) {
            // 32bit modulus
            L_32 = L;
            v1 %= L_32;
            t2_32 = c1 % L_32;
            t1_32 = mulmod_s32(t2_32, a1, L_32);
            t1_32 = mulmod_s32(t1_32, v1, L_32);
            t1_32 = mulmod_s32(b1 - t1_32, v1, L_32) - 2;
            t1_32 = mulmod_s32(t1_32, v1, L_32);
            t1_32 = mulmod_s32(t1_32, t2_32, L_32);
            if (t1_32 < 0) {
                t1_32 += L_32;
            }
            K = t1_32;
        }
        else {
            // 64bit modulus
            K = mulmod_s64(c1, a1, L);
            K = mulmod_s64(K, v1, L);
            K = mulmod_s64(b1 - K, v1, L) - 2;
            K = mulmod_s64(K, v1, L);
            K = mulmod_s64(K, c1, L);
            if (K < 0) {
                K += L;
            }
        }
    }
    else {
        // S = u2 (a SP) + v2 (b^2 - ac)
        temp = (int64_t)SP * (int64_t)a1;
        temp2 = (int64_t)b1 * (int64_t)b1 - (int64_t)a1 * c1;
        if (s64_is_s32(temp) && s64_is_s32(temp2)) {
            // (a*SP) and (b^2-ac) fit into 32bit each
            S = gcdext_s32(&u2, &v2, temp, temp2);
            u2_64 = u2;
        }
        else {
            // 64bit gcdext required
            S = gcdext_s64(&u2_64, &v2_64, temp, temp2);
            v2 = v2_64;
        }

        // N = a/S
        N = a1 / S;

        // L = N a
        L = (int64_t)N * (int64_t)a1;


        // K = -c(v1 u2 a + v2 b) (mod L)
        if (s64_is_s32(L)) {
            // 32bit modulus
            L_32 = L;
            K = u2_64 % L_32;
            K = mulmod_s32(K, v1, L_32);
            K = mulmod_s32(K, a1, L_32);
            t1_32 = mulmod_s32(v2, b1, L_32);
            K = addmod_s32(K, t1_32, L_32);
            t1_32 = (-c1) % L_32;
            K = mulmod_s64(K, t1_32, L_32);
            if (K < 0) {
                K += L_32;
            }
        }
        else {
            // 64bit modulus
            K = mulmod_s64(u2_64, v1, L);
            K = mulmod_s64(K, a1, L);
            temp = mulmod_s64(v2, b1, L);
            K = addmod_s64(K, temp, L);
            K = mulmod_s64(K, -c1, L);
            if (K < 0) {
                K += L;
            }
        }

        // C = Sc
        c1 *= S;
    }

    // Compute NUCOMP termination bound
    B = (int64_t)group->S * (int64_t)a1;
    //B >>= 1;
    //B = sqrt_u64(B);
    B >>= ((msb_u64(B)+1)>>1) + 1; // approximate sqrt, (see above two lines)
    if (B == 0) B = 1;
    

    if (L < B) {
        // compute with regular cubing formula (result will be reduced)

        // T = NK
        T = N * (int32_t)K;

        // R.a = NL
        R->a = N * (int32_t)L;

        // C.b = b + 2 T
        R->b = b1 + (T << 1);

        // C.c = (S c + K (T + b)) / L
        R->c = (K * (T+b1) + c1) / L;
    }
    else {
        // use NUCOMP formulas

        // Execute partial reduction    
        R2_64 = L;
        R1_64 = K;
        gcdext_partial_s64((uint64_t*)&R2_64, (uint64_t*)&R1_64, &C2_64, &C1_64, B);
        R1 = R1_64;
        C1 = C1_64;
        C2 = C2_64;

        // T = N K (mod L)
        T_64 = mulmod_s64(K, N, L);
        if (T_64 < 0) {
            T_64 += L;
        }

        // M1 = (N R1 + T C1) / L
        if (s64_is_s32(T_64)) {
            M1 = muladdmul_s64_4s32(N, R1, T_64, C1) / L;
        }
        else {
            M1 = muladdmuldiv_s64(N, R1, T_64, C1, L);
        }

        // M2 = (R1(b + T) - c S C1) / L
        M2 = muladdmuldiv_s64(T_64+b1, R1, -c1, C1, L);

        // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
        R->a = R1*M1 - C1*M2;
        if (C1 > 0) {
            R->a = -R->a;
        }

        // C.b = 2 (N R1 - C.a C2) / C1 - b (mod 2C.a)
        temp = muladdmul_s64_4s32(N, R1, -R->a, C2);
        temp <<= 1;
        temp /= C1;
        temp -= b1;
        temp2 = (int64_t)R->a<<1;
        temp %= temp2;
        // temp is the closest remainder to zero
        if (R->a < 0) {
            // R->a is negative
            if (temp < R->a) {
                temp -= temp2;
            }
            if (temp > -R->a) {
                temp += temp2;
            }
        }
        else {
            // R->a is positive
            if (temp < -R->a) {
                temp += temp2;
            }
            if (temp > R->a) {
                temp -= temp2;
            }
        }
        R->b = temp;

        // make R->a positive
        if (R->a < 0) {
            R->a = -R->a;
        }

        // C.c = (C.b^2 - D) / 4 C.a
        R->c = s64_qform_c(group, R->a, R->b);
    }

    // normalize and reduce
    s64_qform_reduce(group, R);
}




