/**
 * @file s64_qform.h
 * 
 * Signed 64-bit binary quadratic forms.
 * 
 * These methods are tested to be safe for discriminants up to 59 bits.
 */
#pragma once
#ifndef S64_QFORM__INCLUDED
#define S64_QFORM__INCLUDED

#include <gmp.h>
#include <stdint.h>
#include <stdio.h>

#include "libqform/qform_group.h"

extern const group_cost_t s64_qform_costs;

#define s64_qform_group_max_bits 59

typedef struct {
  int32_t a;
  int32_t b;
  int64_t c;
} s64_qform_t;

typedef struct {
  qform_group_t desc;
  
  int64_t D;
  uint32_t S; // square root of delta
  uint32_t L; // 4th root of delta
} s64_qform_group_t;

void s64_qform_group_init(s64_qform_group_t* group);
static inline void s64_qform_group_clear(s64_qform_group_t* group) {}
void s64_qform_group_set_discriminant(s64_qform_group_t* group, const mpz_t D);
void s64_qform_group_set_discriminant_s64(s64_qform_group_t* group, const int64_t D);

static inline uint32_t s64_qform_hash32(s64_qform_group_t* group, const s64_qform_t* form);
static inline void s64_qform_init(s64_qform_group_t* group, s64_qform_t* R) {}
static inline void s64_qform_clear(s64_qform_group_t* group, s64_qform_t* R) {}
void s64_qform_set_id(s64_qform_group_t* group, s64_qform_t* form);
static inline int s64_qform_is_id(s64_qform_group_t* group, const s64_qform_t* form);
static inline void s64_qform_set(s64_qform_group_t* group, s64_qform_t* R, const s64_qform_t* A);
static inline int s64_qform_equal(s64_qform_group_t* group, const s64_qform_t* A, const s64_qform_t* B);
static inline void s64_qform_inverse(s64_qform_group_t* group, s64_qform_t* form);
void s64_qform_compose(s64_qform_group_t* group, s64_qform_t* R, const s64_qform_t* A, const s64_qform_t* B);
void s64_qform_square(s64_qform_group_t* group, s64_qform_t* R, const s64_qform_t* A);
void s64_qform_cube(s64_qform_group_t* group, s64_qform_t* R, const s64_qform_t* A);
void s64_qform_print(s64_qform_group_t* group, const s64_qform_t* form);

void s64_qform_reduce(s64_qform_group_t* group, s64_qform_t* form);
int s64_qform_is_primeform(s64_qform_group_t* group, s64_qform_t* form, const int p);
static inline int s64_qform_is_ambiguous(s64_qform_group_t* group, const s64_qform_t* form);
int s64_qform_split_ambiguous(s64_qform_group_t* group, mpz_t d, const mpz_t N, const s64_qform_t* form);

/**
 * Inline methods
 */

static inline void s64_qform_set(s64_qform_group_t* group, s64_qform_t* R, const s64_qform_t* A) {
  R->a = A->a;
  R->b = A->b;
  R->c = A->c;
}

static inline void s64_qform_set3(s64_qform_group_t* group, s64_qform_t* R, const int32_t a, const int32_t b, const int64_t c) {
  R->a = a;
  R->b = b;
  R->c = c;
}

static inline int s64_qform_equal(s64_qform_group_t* group, const s64_qform_t* A, const s64_qform_t* B) {
  return A->a == B->a && A->b == B->b;  // c is computed from a and b
}

static inline uint32_t s64_qform_hash32(s64_qform_group_t* group, const s64_qform_t* form) {
  // Magic number is smallest prime larger than (1+sqrt(5))/2 * 2^31,
  // which is the golden ratio.
  // NOTE: It's important that we only use the a and c components, since
  //       b is redundant, but more importantly, this means that
  //       the form and its inverse hash to the same value
  //       and will get picked up by collision detection routines.
  const uint32_t magic = 3474701543UL;
  return (((uint32_t)form->a * magic) + (uint32_t)form->c) * magic;
}

static inline int s64_qform_is_id(s64_qform_group_t* group, const s64_qform_t* form) {
  return form->a == 1;
}

static inline void s64_qform_inverse(s64_qform_group_t* group, s64_qform_t* form) {
  if (form->a != form->b && form->a != form->c) {
    form->b = -form->b;
  }
}

static inline int s64_qform_is_ambiguous(s64_qform_group_t* group, const s64_qform_t* form) {
  return form->a > 1 && (form->b == 0 || form->a == form->b || form->c == form->a);
}

#endif  // S64_QFORM__INCLUDED

