/**
 * @file qform_group.h
 * 
 * Binary Quadratic Form Group Descriptor.
 */
#pragma once
#ifndef QFORM_GROUP__INCLUDED
#define QFORM_GROUP__INCLUDED

#include <gmp.h>
#include <stdlib.h>

#include "liboptarith/closest_23.h"
#include "liboptarith/group.h"

typedef void qform_t;
typedef struct qform_group_struct qform_group_t;

typedef void qform_group_clear_f(qform_group_t* group);
typedef void qform_group_set_discriminant_f(qform_group_t* group, const mpz_t D);
typedef void qform_reduce_f(qform_group_t* group, qform_t* form);
typedef int qform_is_primeform_f(qform_group_t* group, qform_t* form, const int p);
typedef int qform_is_ambiguous_f(qform_group_t* group, const qform_t* form);
typedef int qform_split_ambiguous_f(qform_group_t* group, mpz_t d, const mpz_t N, const qform_t* form);

struct qform_group_struct {
  // 'group' must be the first member of the struct.
  group_t group;
  
  int discriminant_max_bits;
  qform_group_clear_f* clear;
  qform_group_set_discriminant_f* set_discriminant;
  qform_reduce_f* reduce;
  qform_is_primeform_f* is_primeform;
  qform_is_ambiguous_f* is_ambiguous;
  qform_split_ambiguous_f* split_ambiguous;

  // Precomputed 2,3 representations for uint16_t.
  int* pow_rep_sizes;
  factored_two_three_term16_t** pow_reps;
};

/**
 * Compute a random prime form from the prime list.
 */
void qform_random_primeform(qform_group_t* group, qform_t* form);

/**
 * Compute the next prime form and return the index of the prime.
 * The return value is guaranteed to be greater than 'prime_index'.
 */
int qform_next_primeform(qform_group_t* group,
                         qform_t* form,
                         int prime_index);

#endif 


