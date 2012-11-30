/**
 * GNU MPZ version of qforms
 */
#include <limits.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>

#include "liboptarith/math_mpz.h"
#include "liboptarith/mpz_xgcd.h"
#include "liboptarith/primes.h"
#include "liboptarith/sqrtmodp_list.h"
#include "libqform/mpz_qform.h"

const group_cost_t mpz_qform_costs = {
  1.187,
  1,
  1.624
};

/**
 * Initialize a qform group structure
 */
void mpz_qform_group_init(mpz_qform_group_t* group) {
  const int nbits = 256;
  
  // initialize group description
  group->desc.group.elem_init = (group_elem_init_f*)&mpz_qform_init;
  group->desc.group.elem_clear = (group_elem_clear_f*)&mpz_qform_clear;
  group->desc.group.elem_size = sizeof(mpz_qform_t);
  group->desc.group.hash32 = (group_hash32_f*)&mpz_qform_hash32;
  group->desc.group.set_id = (group_set_id_f*)&mpz_qform_set_id;
  group->desc.group.is_id = (group_is_id_f*)&mpz_qform_is_id;
  group->desc.group.set = (group_set_f*)&mpz_qform_set;
  group->desc.group.equal = (group_equal_f*)&mpz_qform_equal;
  group->desc.group.inverse = (group_inverse_f*)&mpz_qform_inverse;
  group->desc.group.compose = (group_compose_f*)&mpz_qform_compose;
  group->desc.group.square = (group_square_f*)&mpz_qform_square;
  group->desc.group.cube = (group_cube_f*)&mpz_qform_cube;
  group->desc.group.print = (group_print_f*)&mpz_qform_print;
  
  group->desc.discriminant_max_bits = INT_MAX;
  group->desc.clear = (qform_group_clear_f*)&mpz_qform_group_clear;
  group->desc.set_discriminant = (qform_group_set_discriminant_f*)&mpz_qform_group_set_discriminant;
  group->desc.reduce = (qform_reduce_f*)&mpz_qform_reduce;
  group->desc.is_primeform = (qform_is_primeform_f*)&mpz_qform_is_primeform;
  group->desc.is_ambiguous = (qform_is_ambiguous_f*)&mpz_qform_is_ambiguous;
  group->desc.split_ambiguous = (qform_split_ambiguous_f*)&mpz_qform_split_ambiguous;
  
  mpz_init2(group->D, nbits);
  mpz_init2(group->S, nbits);
  mpz_init2(group->L, nbits);
  mpz_init2(group->t, nbits);
  mpz_init2(group->m, nbits);
  mpz_init2(group->n, nbits);
  
  // temporaries for reduction
  mpz_init2(group->reduce.x, nbits);
  mpz_init2(group->reduce.q, nbits);
  mpz_init2(group->reduce.r, nbits);
  
  // temporaries for composition (nucomp, nudupl, nucube)
  mpz_init2(group->compose.a1, nbits);
  mpz_init2(group->compose.b1, nbits);
  mpz_init2(group->compose.c1, nbits);
  mpz_init2(group->compose.a2, nbits);
  mpz_init2(group->compose.b2, nbits);
  mpz_init2(group->compose.c2, nbits);
  mpz_init2(group->compose.ss, nbits);
  mpz_init2(group->compose.m, nbits);
  mpz_init2(group->compose.SP, nbits);
  mpz_init2(group->compose.S, nbits);
  mpz_init2(group->compose.v1, nbits);
  mpz_init2(group->compose.u2, nbits);
  mpz_init2(group->compose.v2, nbits);
  mpz_init2(group->compose.K, nbits);
  mpz_init2(group->compose.T, nbits);
  mpz_init2(group->compose.temp, nbits);
  mpz_init2(group->compose.temp2, nbits);
  mpz_init2(group->compose.R1, nbits);
  mpz_init2(group->compose.R2, nbits);
  mpz_init2(group->compose.C1, nbits);
  mpz_init2(group->compose.C2, nbits);
  mpz_init2(group->compose.M1, nbits);
  mpz_init2(group->compose.M2, nbits);
  mpz_init2(group->compose.N, nbits);
  mpz_init2(group->compose.L, nbits);
  mpz_init2(group->compose.B, nbits);
  mpz_xgcd_init(&group->compose.gcd, nbits);     
}

void mpz_qform_group_clear(mpz_qform_group_t* group) {
  mpz_clear(group->D);
  mpz_clear(group->S);
  mpz_clear(group->L);
  mpz_clear(group->t);
  mpz_clear(group->m);
  mpz_clear(group->n);
  
  // temporaries for reduction
  mpz_clear(group->reduce.x);
  mpz_clear(group->reduce.q);
  mpz_clear(group->reduce.r);
  
  // temporaries for composition (nucomp, nudupl, nucube)
  mpz_clear(group->compose.a1);
  mpz_clear(group->compose.b1);
  mpz_clear(group->compose.c1);
  mpz_clear(group->compose.a2);
  mpz_clear(group->compose.b2);
  mpz_clear(group->compose.c2);
  mpz_clear(group->compose.ss);
  mpz_clear(group->compose.m);
  mpz_clear(group->compose.SP);
  mpz_clear(group->compose.S);
  mpz_clear(group->compose.v1);
  mpz_clear(group->compose.u2);
  mpz_clear(group->compose.v2);
  mpz_clear(group->compose.K);
  mpz_clear(group->compose.T);
  mpz_clear(group->compose.temp);
  mpz_clear(group->compose.temp2);
  mpz_clear(group->compose.R1);
  mpz_clear(group->compose.R2);
  mpz_clear(group->compose.C1);
  mpz_clear(group->compose.C2);
  mpz_clear(group->compose.M1);
  mpz_clear(group->compose.M2);
  mpz_clear(group->compose.N);
  mpz_clear(group->compose.L);
  mpz_clear(group->compose.B);
  mpz_xgcd_clear(&group->compose.gcd);     
}

void mpz_qform_group_set_discriminant(mpz_qform_group_t* group, const mpz_t D) {
  mpz_set(group->D, D);
  mpz_abs(group->S, D);  // Positive sqrt, since negative is undefined.
  mpz_sqrt(group->S, group->S);
  mpz_sqrt(group->L, group->S);
}

// compute c = (b^2-D)/(4a)
void mpz_qform_c(mpz_qform_group_t* group, mpz_t c, mpz_t a, mpz_t b) {
  mpz_mul (c, b, b);
  mpz_sub (c, c, group->D);
  mpz_fdiv_q_2exp (c, c, 2);
  mpz_divexact (c, c, a);
}

void mpz_qform_set_id(mpz_qform_group_t* group, mpz_qform_t* form) {
  mpz_set_ui(form->a, 1);
  
  // D (mod 4)
  if ((group->D->_mp_d[0] & 3) == 0) {
    mpz_set_ui(form->b, 0);
  } else {
    mpz_set_ui(form->b, 1);
  }
  
  mpz_qform_c(group, form->c, form->a, form->b);
}

/**
 * If N is split by the form (a,b,C), then d is a factor and 1 is returned
 * otherwise 0 is returned.
 * Assumes that form is an ambiguous form.
 */
int mpz_qform_split_ambiguous(mpz_qform_group_t* group, mpz_t out_d, const mpz_t in_N, const mpz_qform_t* form) {
  // ambiguous forms are of three types
  // (a, 0, c), (a, a, c), and (a, b, a)
  if (mpz_cmp_ui(form->b, 0) == 0 || mpz_cmp(form->a, form->b) == 0) {
    // first two forms (a, 0, c) and (a, a, c)
    // (a, 0, c) => -4ac = D => a | D
    // (a, a, c) => a^2-4ac = a(a-4c) = D => a | d
    mpz_set(group->m, form->a);
  } else {
    // final form (a, b, a)
    // (a, b, a) => b^2 - 4a^2 = (b+2a)(b-2a) = D
    // pick m = min |b +/- 2a|
    mpz_mul_2exp(group->m, form->a, 1);
    if (form->b->_mp_size < 0) {
      mpz_add(group->m, group->m, form->b);
    } else {
      mpz_sub(group->m, form->b, group->m);
    }
  }
  
  mpz_abs(group->m, group->m);
  mpz_gcd(out_d, group->m, in_N);

  // return true if this is a non-trivial root
  return (mpz_cmp_ui(out_d, 1) > 0 && mpz_cmp(out_d, in_N) < 0);
}

/**
 * Check for a primeform where a=p.
 * Tests b = +/- sqrt(D) mod p.
 * Note: We let -0 mod p = p.  This way we try ambiguous forms.
 * TODO: sqrtmodp should work for any inputs
 */
int mpz_qform_is_primeform(mpz_qform_group_t* group, mpz_qform_t* form, int p) {
  int Dmodp;
  const short* sqrtp;
  
  if (p > sqrtmodp_maxp)
    return 0; // p is too large for the table
  
  sqrtp = sqrtmodp[p];
  if (sqrtp == 0)
    return 0; // p was not prime
    
  // a = p
  mpz_set_s64(form->a, p);
  
  // b = sqrt(D) (mod p)
  Dmodp = mpz_fdiv_ui(group->D, p);
  if (Dmodp < 0) {
    Dmodp += p;
  }
  if (Dmodp == 0) {
    // p divides the discriminant
    return 0;
  }
  if (sqrtp[Dmodp] == -1) {
    return 0; // D does not have a sqrt (mod p)
  }
  mpz_set_s64(form->b, sqrtp[Dmodp]);
  
  // We know that p | b^2-D for +/- b.
  // Compute c = (b^2 - d) / (4a) if possible.
  if (p == 2) {
    // special case if p == 2
    mpz_mul(form->c, form->b, form->b);
    mpz_sub(form->c, form->c, group->D);

    if ((form->c->_mp_d[0] & 7) == 0) {
      // 4a | b^2-D, a=2, so we divide by 8
      mpz_fdiv_q_2exp(form->c, form->c, 3);
      return 1;
    }
    if (sqrtp[Dmodp] == 1) // form->b == 1
      return 0;
    
    // try -b = 2
    mpz_set_s64(form->b, 2);
    mpz_set_s64(form->c, 4);
    mpz_sub(form->c, form->c, group->D);
    if ((form->c->_mp_d[0] & 7) != 0)
      return 0;
    
    // 4a | b^2-D, a=2, so we divide by 8
    mpz_fdiv_q_2exp(form->c, form->c, 3);
    return 1;
  }

  // p != 2
  mpz_mul(form->c, form->b, form->b);
  mpz_sub(form->c, form->c, group->D);
  
  if ((form->c->_mp_d[0] & 3) == 0) {
    // 4a | b^2-D
    // divide by 4
    mpz_fdiv_q_2exp(form->c, form->c, 2);
    
    // divide by a
    mpz_divexact_ui(form->c, form->c, p);
    return 1;
  }

  // try -b
  mpz_set_s64(form->b, p-sqrtp[Dmodp]);
  
  // make sure that 4 | b^2-D
  mpz_mul(form->c, form->b, form->b);
  mpz_sub(form->c, form->c, group->D);
  if ((form->c->_mp_d[0] & 3) != 0) {
    return 0;
  }
  
  // divide by 4
  mpz_fdiv_q_2exp(form->c, form->c, 2);
  
  // divide by a
  mpz_divexact_ui(form->c, form->c, p);
  
  return 1;
}

/**
 * Form is assumed to be positive definite, primitive, and the discriminant is assumed to be negative
 * Algorithm 5.4.2 in Cohen  "A Course in Computational Algebraic Number Theory" p. 243
 */
void mpz_qform_reduce(mpz_qform_group_t* group, mpz_qform_t* form) {
  register mpz_qform_reduce_t* reduce = &group->reduce;
  int i = 0;
  
  // First make sure form is positive definite
  while (1) {
    if (mpz_cmp(form->a, form->c) > 0) {
      mpz_set(reduce->x, form->a);
      mpz_set(form->a, form->c);
      mpz_set(form->c, reduce->x);
      mpz_neg(form->b, form->b);
    }
    i = mpz_cmpabs(form->b, form->a);
    if (i > 0 || (i == 0 && mpz_sgn(form->b) == -1)) {
      mpz_mul_2exp(reduce->x, form->a, 1);
      mpz_fdiv_qr(reduce->q, reduce->r, form->b, reduce->x);
      if (mpz_cmp(reduce->r, form->a) > 0) {
	mpz_sub(reduce->r, reduce->r, reduce->x);
	mpz_add_ui(reduce->q, reduce->q, 1);
      }
      mpz_add(reduce->x, form->b, reduce->r);
      mpz_mul(reduce->x, reduce->x, reduce->q);
      mpz_fdiv_q_2exp(reduce->x, reduce->x, 1);
      mpz_sub(form->c, form->c, reduce->x);
      mpz_set(form->b, reduce->r);
    } else {
      break;
    }
  }
  if (mpz_cmp(form->a, form->c) == 0 && mpz_sgn(form->b) < 0) {
    mpz_neg(form->b, form->b);
  }
}

/**
 * NUCOMP algorithm. Adapted from "Solving the Pell Equation"
 * by Michael J. Jacobson, Jr. and Hugh C. Williams.
 * http://www.springer.com/mathematics/numbers/book/978-0-387-84922-5
 */
void mpz_qform_compose(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A, const mpz_qform_t* B) {
  register mpz_qform_compose_t* comp = &group->compose;
  
  if (mpz_cmp(A->a, B->a) > 0) {
    mpz_set(comp->a1, A->a);
    mpz_set(comp->b1, A->b);
    mpz_set(comp->a2, B->a);
    mpz_set(comp->b2, B->b);
    mpz_set(comp->c2, B->c);
  } else {
    mpz_set(comp->a1, B->a);
    mpz_set(comp->b1, B->b);
    mpz_set(comp->a2, A->a);
    mpz_set(comp->b2, A->b);
    mpz_set(comp->c2, A->c);
  }
  
  // s = (b1 + b2)/2, m = (b1 - b2)/2
  mpz_add(comp->ss, comp->b1, comp->b2);
  mpz_fdiv_q_2exp(comp->ss, comp->ss, 1);
  
  mpz_sub(comp->m, comp->b1, comp->b2);
  mpz_fdiv_q_2exp(comp->m, comp->m, 1);
  
  // solve SP = v1 a2 + u1 a1 (only need v1)
  mpz_gcdext(comp->SP, comp->v1, 0, comp->a2, comp->a1);
  
  // K = v1 (b1 - b2) / 2 (mod L)
  mpz_mul(comp->K, comp->m, comp->v1);
  mpz_fdiv_r(comp->K, comp->K, comp->a1);
  
  if (mpz_cmp_ui(comp->SP, 1) != 0) {
    mpz_gcdext(comp->S, comp->u2, comp->v2, comp->SP, comp->ss);
    
    // K = u2 K - v2 c2 (mod L)
    mpz_mul(comp->K, comp->K, comp->u2);
    mpz_mul(comp->temp, comp->v2, comp->c2);
    mpz_sub(comp->K, comp->K, comp->temp);
    
    if (mpz_cmp_ui(comp->S, 1) != 0) {
      mpz_divexact(comp->a1, comp->a1, comp->S);
      mpz_divexact(comp->a2, comp->a2, comp->S);
      mpz_mul(comp->c2, comp->c2, comp->S);
    }
    
    mpz_fdiv_r(comp->K, comp->K, comp->a1);
  }
  
  // N = a2;  L = a1;
  
  // check if NUCOMP steps are required
  if (mpz_cmp(comp->a1, group->L) < 0) {
    // compute with regular multiplication formula (result will be reduced)
    
    // T = NK
    mpz_mul(comp->T, comp->a2, comp->K);
    
    // C.a = A.a B.a / d^2 = NL
    mpz_mul(R->a, comp->a2, comp->a1);
    
    // C.b = b2 + 2 a2 K = b2 + 2 T
    mpz_mul_2exp(R->b, comp->T, 1);
    mpz_add(R->b, R->b, comp->b2);
    
    // C.c = (S c2 + K (b2 + T)) / L;
    mpz_add(R->c, comp->b2, comp->T);
    mpz_mul(R->c, R->c, comp->K);
    mpz_add(R->c, R->c, comp->c2);
    mpz_divexact(R->c, R->c, comp->a1);
    
    mpz_qform_reduce(group, R);
  } else {
    // use NUCOMP formulas
    
    // Execute partial reduction
    mpz_set(comp->R2, comp->a1);
    mpz_set(comp->R1, comp->K);
    mpz_xgcd_partial(&comp->gcd, comp->R2, comp->R1, comp->C2, comp->C1, group->L);
    
    // M1 = (N R1 + (b1 - b2) C1 / 2) / L  (T = N R1)
    mpz_mul(comp->T, comp->a2, comp->R1);
    mpz_mul(comp->M1, comp->m, comp->C1);
    mpz_add(comp->M1, comp->M1, comp->T);
    mpz_divexact(comp->M1, comp->M1, comp->a1);
    
    // M2 = (R1(b1 + b2)/2 - c2 S C1) / L
    mpz_mul(comp->M2, comp->ss, comp->R1);
    mpz_mul(comp->temp, comp->c2, comp->C1);
    mpz_sub(comp->M2, comp->M2, comp->temp);
    mpz_divexact(comp->M2, comp->M2, comp->a1);
    
    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    mpz_mul(R->a, comp->R1, comp->M1);
    mpz_mul(comp->temp, comp->C1, comp->M2);
    if (mpz_sgn(comp->C1) < 0) {
      mpz_sub(R->a, R->a, comp->temp);
    } else {
      mpz_sub(R->a, comp->temp, R->a);
    }
    
    // C.b = 2 (N R1 - C.a C2) / C1 - b2 (mod 2a)
    mpz_mul(R->b, R->a, comp->C2);
    mpz_sub(R->b, comp->T, R->b);
    mpz_mul_2exp(R->b, R->b, 1);
    mpz_divexact(R->b, R->b, comp->C1);
    mpz_sub(R->b, R->b, comp->b2);
    mpz_mul_2exp(comp->temp, R->a, 1);
    mpz_fdiv_r(R->b, R->b, comp->temp);
    
    // C.c = (C.b^2 - Delta) / 4 C.a
    mpz_mul(R->c, R->b, R->b);
    mpz_sub(R->c, R->c, group->D);
    mpz_fdiv_q_2exp(R->c, R->c, 2);
    mpz_divexact(R->c, R->c, R->a);
    
    if (mpz_sgn(R->a) < 0) {
      mpz_neg(R->a, R->a);
      mpz_neg(R->c, R->c);
    }
    
    mpz_qform_reduce(group, R);
  }
}

/**
 * NUDUPL. Simplified from compose above.
 */
void mpz_qform_square(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A) {
  register mpz_qform_compose_t* comp = &group->compose;
  
  mpz_set(comp->a1, A->a);
  mpz_set(comp->b1, A->b);
  mpz_set(comp->c1, A->c);
  
  // solve S = v1 b1 + u1 a1 (only need v1)
  mpz_gcdext(comp->S, comp->v1, 0, comp->b1, comp->a1);
  
  // K = -v1 c1 (mod L)
  mpz_mul(comp->K, comp->v1, comp->c1);
  mpz_neg(comp->K, comp->K);
  
  if (mpz_cmp_ui(comp->S, 1) != 0) {
    mpz_divexact(comp->a1, comp->a1, comp->S);
    mpz_mul(comp->c1, comp->c1, comp->S);
  }
  
  mpz_fdiv_r(comp->K, comp->K, comp->a1);
  
  // N = L = a1;
  
  // check if NUCOMP steps are required
  if (mpz_cmp(comp->a1, group->L) < 0) {
    // compute with regular squaring formula (result will be reduced)
    
    // T = NK
    mpz_mul(comp->T, comp->a1, comp->K);
    
    // C.a = A.a^2 / S^2 = N^2
    mpz_mul(R->a, comp->a1, comp->a1);
    
    // C.b = b1 + 2 a1 K = b1 + 2 T
    mpz_mul_2exp(R->b, comp->T, 1);
    mpz_add(R->b, R->b, comp->b1);
    
    // C.c = (S c1 + K (b1 + T)) / L;
    mpz_add(R->c, comp->b1, comp->T);
    mpz_mul(R->c, R->c, comp->K);
    mpz_add(R->c, R->c, comp->c1);
    mpz_divexact(R->c, R->c, comp->a1);
    
    mpz_qform_reduce(group, R);
  } else {
    // use NUCOMP formulas
    
    // Execute partial reduction
    mpz_set(comp->R2, comp->a1);
    mpz_set(comp->R1, comp->K);
    mpz_xgcd_partial(&comp->gcd, comp->R2, comp->R1, comp->C2, comp->C1, group->L);
    
    // M1 = R1
    
    // M2 = (R1 b1 - c1 S C1) / L
    mpz_mul(comp->M2, comp->R1, comp->b1);
    mpz_mul(comp->temp, comp->c1, comp->C1);
    mpz_sub(comp->M2, comp->M2, comp->temp);
    mpz_divexact(comp->M2, comp->M2, comp->a1);
    
    // C.a = (-1)^(i-1) (R1^2 - C1 M2)
    mpz_mul(R->a, comp->R1, comp->R1);
    mpz_mul(comp->temp, comp->C1, comp->M2);
    if (mpz_sgn(comp->C1) < 0) {
      mpz_sub(R->a, R->a, comp->temp);
    } else {
      mpz_sub(R->a, comp->temp, R->a);
    }
    
    // C.b = 2 (N R1 - C.a C2) / C1 - b1 (mod 2a)
    mpz_mul(R->b, R->a, comp->C2);
    mpz_mul(comp->temp, comp->a1, comp->R1);
    mpz_sub(R->b, comp->temp, R->b);
    mpz_mul_2exp(R->b, R->b, 1);
    mpz_divexact(R->b, R->b, comp->C1);
    mpz_sub(R->b, R->b, comp->b1);
    mpz_mul_2exp(comp->temp, R->a, 1);
    mpz_fdiv_r(R->b, R->b, comp->temp);
    
    // C.c = (C.b^2 - Delta) / 4 C.a
    mpz_mul(R->c, R->b, R->b);
    mpz_sub(R->c, R->c, group->D);
    mpz_fdiv_q_2exp(R->c, R->c, 2);
    mpz_divexact(R->c, R->c, R->a);
    
    if (mpz_sgn(R->a) < 0) {
      mpz_neg(R->a, R->a);
      mpz_neg(R->c, R->c);
    }
    
    mpz_qform_reduce(group, R);
  }
}

/**
 * Computes a reduced ideal equivalent to the cube of an ideal.
 * Adapted from "Fast Ideal Cubing in Imaginary Quadratic Number
 * and Function Fields" by Laurent Imbert, Michael J. Jacobson, Jr. and
 * Arthur Schmidt.
 * www.lirmm.fr/~imbert/pdfs/cubing_amc_2010.pdf
 */
void mpz_qform_cube(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A) {
  register mpz_qform_compose_t* comp = &group->compose;
  
  // special case when (a,b) is an ambiguous form
  if (mpz_qform_is_ambiguous(group, A)) {
    // cubing an ambiguous form gets us the form back
    mpz_qform_set(group, R, A);
    mpz_qform_reduce(group, R);
    return;
  }
  
  mpz_set(comp->b1, A->b);
  mpz_set(comp->c1, A->c);
  
  // solve SP = v1 b + u1 a (only need v1)
  mpz_gcdext(comp->SP, comp->v1, 0, A->b, A->a);
  
  if (mpz_cmp_ui(comp->SP, 1) == 0) {
    // N = a
    mpz_set(comp->N, A->a);
    
    // L = a^2
    mpz_mul(comp->L, A->a, A->a);
    
    // K = c v1 (v1(b - a c v1) - 2) mod L
    mpz_fdiv_r(comp->v1, comp->v1, comp->L);
    mpz_mul(comp->temp, comp->v1, comp->c1);
    mpz_mul(comp->K, comp->temp, A->a);
    mpz_sub(comp->K, A->b, comp->K);
    mpz_mul(comp->K, comp->K, comp->v1);
    mpz_sub_ui(comp->K, comp->K, 2);
    mpz_mulm(comp->K, comp->K, comp->temp, comp->L);
  } else {
    // S = u2 (a SP) + v2 (b^2 - ac)
    mpz_mul(comp->SP, comp->SP, A->a);
    
    mpz_mul(comp->temp, A->b, A->b);
    mpz_mul(comp->T, A->a, comp->c1);
    mpz_sub(comp->temp, comp->temp, comp->T);
    
    mpz_gcdext(comp->S, comp->u2, comp->v2, comp->SP, comp->temp);
    
    // N = a/S
    mpz_divexact(comp->N, A->a, comp->S);
    
    // L = N a
    mpz_mul(comp->L, comp->N, A->a);
    
    // K = -c(v1 u2 a + v2 b) mod L
    mpz_fdiv_r(comp->v2, comp->v2, comp->L);
    mpz_mul(comp->K, comp->v1, comp->u2);
    mpz_mul(comp->K, comp->K, A->a);
    mpz_mul(comp->temp, comp->v2, A->b);
    mpz_add(comp->K, comp->K, comp->temp);
    mpz_mul(comp->K, comp->K, comp->c1);
    mpz_sub(comp->K, comp->L, comp->K);
    mpz_fdiv_r(comp->K, comp->K, comp->L);
    
    // C = Sc
    mpz_mul(comp->c1, comp->c1, comp->S);
  }

  // Compute NUCOMP termination bound
  mpz_mul(comp->B, group->S, A->a);
  // use a sqrt approximation
  // (TODO: is approximation what we want?)
  mpz_fdiv_q_2exp(comp->B, comp->B, (mpz_sizeinbase(comp->B, 2)>>1)+1);
  //mpz_fdiv_q_2exp(comp->B, comp->B, 1);
  //mpz_sqrt(comp->B, comp->B);
  
  if (mpz_cmp(comp->L, comp->B) < 0) {
    // compute with regular cubing formula (result will be reduced)
    
    // T = NK
    mpz_mul(comp->T, comp->N, comp->K);
    
    // C.a = N L
    mpz_mul(R->a, comp->N, comp->L);
    
    // C.b = b + 2 T
    mpz_mul_2exp(R->b, comp->T, 1);
    mpz_add(R->b, R->b, comp->b1);
    
    // C.c = (S c + K (T + b)) / L
    mpz_add(R->c, comp->T, comp->b1);
    mpz_mul(R->c, R->c, comp->K);
    mpz_add(R->c, R->c, comp->c1);
    mpz_divexact(R->c, R->c, comp->L);
    
    mpz_qform_reduce(group, R);
  } else {
    // use NUCOMP formulas
    
    // Execute partial reduction
    mpz_set(comp->R2, comp->L);
    mpz_set(comp->R1, comp->K);
    mpz_xgcd_partial(&comp->gcd, comp->R2, comp->R1, comp->C2, comp->C1, comp->B);
    
    // T = N K
    mpz_mulm(comp->T, comp->N, comp->K, comp->L);
    
    // M1 = (N R1 + T C1) / L  (temp = N R1)
    mpz_mul(comp->temp, comp->N, comp->R1);
    mpz_mul(comp->M1, comp->T, comp->C1);
    mpz_add(comp->M1, comp->M1, comp->temp);
    mpz_divexact(comp->M1, comp->M1, comp->L);
    
    // M2 = (R1(b + T) - c S C1) / L
    mpz_add(comp->M2, A->b, comp->T);
    mpz_mul(comp->M2, comp->M2, comp->R1);
    mpz_mul(comp->temp2, comp->c1, comp->C1);
    mpz_sub(comp->M2, comp->M2, comp->temp2);
    mpz_divexact(comp->M2, comp->M2, comp->L);
    
    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    mpz_mul(R->a, comp->R1, comp->M1);
    mpz_mul(comp->temp2, comp->C1, comp->M2);
    if (mpz_sgn(comp->C1) < 0) {
      mpz_sub(R->a, R->a, comp->temp2);
    } else {
      mpz_sub(R->a, comp->temp2, R->a);
    }
    
    // C.b = 2 (N R1 - C.a C2) / C1 - b (mod 2C.a)
    mpz_mul(R->b, R->a, comp->C2);
    mpz_sub(R->b, comp->temp, R->b);
    mpz_mul_2exp(R->b, R->b, 1);
    mpz_divexact(R->b, R->b, comp->C1);
    mpz_sub(R->b, R->b, comp->b1);
    mpz_mul_2exp(comp->temp, R->a, 1);
    mpz_fdiv_r(R->b, R->b, comp->temp);
    
    // C.c = (C.b^2 - Delta) / 4 C.a
    mpz_mul(R->c, R->b, R->b);
    mpz_sub(R->c, R->c, group->D);
    mpz_fdiv_q_2exp(R->c, R->c, 2);
    mpz_divexact(R->c, R->c, R->a);
    
    if (mpz_sgn(R->a) < 0) {
      mpz_neg(R->a, R->a);
      mpz_neg(R->c, R->c);
    }
    
    // normalize and reduce
    mpz_qform_reduce(group, R);
  }
}


