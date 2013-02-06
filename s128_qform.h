#pragma once
#ifndef S128_QFORM__INCLUDED
#define S128_QFORM__INCLUDED

#include <gmp.h>
#include <stdint.h>

#include "liboptarith/s128_t.h"
#include "libqform/qform_group.h"

extern const group_cost_t s128_qform_costs;

#define s128_qform_group_max_bits 118

typedef struct {
  int64_t a;
  int64_t b;
  s128_t c;
} s128_qform_t;

typedef struct {
  // must be the first member of the struct
  qform_group_t desc;
  
  s128_t D;
  uint64_t S; // square root of delta
  uint64_t L; // 4th root of delta
  
  mpz_t tmp;  // used to promote intermediates to avoid overflow
  mpz_t tmp2;
  mpz_t tmp3;
} s128_qform_group_t;

void s128_qform_group_init(s128_qform_group_t* group);
void s128_qform_group_clear(s128_qform_group_t* group);
void s128_qform_group_set_discriminant(s128_qform_group_t* group, const mpz_t D);
void s128_qform_group_set_discriminant_s128(s128_qform_group_t* group, const s128_t* D);

static inline void s128_qform_init(s128_qform_group_t* group, s128_qform_t* R) {}
static inline void s128_qform_clear(s128_qform_group_t* group, s128_qform_t* R) {}
static inline uint32_t s128_qform_hash32(s128_qform_group_t* group, const s128_qform_t* form);
void s128_qform_set_id(s128_qform_group_t* group, s128_qform_t* form);
static inline int s128_qform_is_id(s128_qform_group_t* group, const s128_qform_t* form);
static inline void s128_qform_set(s128_qform_group_t* group, s128_qform_t* R, const s128_qform_t* A);
static inline int s128_qform_equal(s128_qform_group_t* group, const s128_qform_t* A, const s128_qform_t* B);
static inline void s128_qform_inverse(s128_qform_group_t* group, s128_qform_t* form);
void s128_qform_compose(s128_qform_group_t* group, s128_qform_t* R, const s128_qform_t* A, const s128_qform_t* B);
void s128_qform_square(s128_qform_group_t* group, s128_qform_t* R, const s128_qform_t* A);
void s128_qform_cube(s128_qform_group_t* group, s128_qform_t* R, const s128_qform_t* A);
void s128_qform_print(s128_qform_group_t* group, const s128_qform_t* form);

void s128_qform_reduce(s128_qform_group_t* group, s128_qform_t* form);
int s128_qform_is_primeform(s128_qform_group_t* group, s128_qform_t* form, const int p);
static inline int s128_qform_is_ambiguous(s128_qform_group_t* group, const s128_qform_t* form);
int s128_qform_split_ambiguous(s128_qform_group_t* group, mpz_t d, const mpz_t N, const s128_qform_t* form);

/**
 * Inline methods
 */
static inline void s128_qform_set(s128_qform_group_t* group, s128_qform_t* R, const s128_qform_t* A) {
  R->a = A->a;
  R->b = A->b;
  R->c = A->c;
}

static inline int s128_qform_equal(s128_qform_group_t* group, const s128_qform_t* A, const s128_qform_t* B) {
  return A->a == B->a && A->b == B->b && is_equal_s128_s128(&A->c, &B->c);
}

static inline void s128_qform_set3(s128_qform_group_t* group, s128_qform_t* R, const int64_t a, const int64_t b, const s128_t* c) {
  R->a = a;
  R->b = b;
  R->c = *c;
}

static inline uint32_t s128_qform_hash32(s128_qform_group_t* group, const s128_qform_t* form) {
  // Magic number is largest 32-bit unsigned prime.
  // It's important that we only use the a and c components, since
  // b is redundant, but more importantly, this means that
  // the form and its inverse hash to the same value
  // and will get picked up by collision detection routines.
  const uint32_t magic = 4294967291UL;
  return (((uint32_t)form->a * magic) + (uint32_t)form->c.v0) * magic;
}

static inline int s128_qform_is_id(s128_qform_group_t* group, const s128_qform_t* form) {
  return form->a == 1;
}

static inline void s128_qform_inverse(s128_qform_group_t* group, s128_qform_t* form) {
  if (form->a != form->b && cmp_s128_s64(&form->c, form->a) != 0) {
    form->b = -form->b;
  }
}

static inline int s128_qform_is_ambiguous(s128_qform_group_t* group, const s128_qform_t* form) {
  return form->a > 1 && (form->b == 0 || form->a == form->b || cmp_s64_s128(form->a, &form->c) == 0);
}

#endif  // S128_QFORM__INCLUDED

