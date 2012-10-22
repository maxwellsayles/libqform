/**
 * Generic qform class.
 * Dynamically chooses the best set of functions
 * depending on the size of the discriminant.
 *
 * We want to do allocation and initialization once only.
 * Therefore, the generic qform code should be more like the sanity qform code
 * in that it allocates and initializes each of s64, s128, and mpz
 * and then dispatches to the desired routine based on the size of the
 * discriminant.  There is a little bit of overhead in terms of memory and
 * speed since each qform must allocate memory for s64, s128, and mpz types
 * and the code for each operation must branch depending on the size of the
 * discriminant.
 */

#pragma once
#ifndef GEN_QFORM_INCLUDED
#define	GEN_QFORM_INCLUDED

#include <inttypes.h>
#include <limits.h>

#include "dbreps/mpz_pow_reps.h"
#include "dbreps/s64_pow_reps.h"
#include "dbreps/s128_pow_reps.h"
#include "liboptarith/group_pow.h"
#include "libqform/mpz_qform.h"
#include "libqform/qform_group.h"
#include "libqform/s128_qform.h"
#include "libqform/s64_qform.h"

typedef struct {
  s64_qform_t s64_qform;
  s128_qform_t s128_qform;
  mpz_qform_t mpz_qform;
} gen_qform_t;

typedef struct {
  qform_group_t desc;
  
  s64_qform_group_t s64_qform_group;
  s128_qform_group_t s128_qform_group;
  mpz_qform_group_t mpz_qform_group;
  
  int logD;
} gen_qform_group_t;
  
static inline void gen_qform_group_init(gen_qform_group_t* this);
static inline void gen_qform_group_clear(gen_qform_group_t* this);

static inline void gen_qform_init(gen_qform_group_t* group, gen_qform_t* this);
static inline void gen_qform_clear(gen_qform_group_t* group, gen_qform_t* this);
static inline uint32_t gen_qform_hash32(gen_qform_group_t* group, const gen_qform_t* form);
static inline void gen_qform_set_id(gen_qform_group_t* group, gen_qform_t* form);
static inline int gen_qform_is_id(gen_qform_group_t* group, const gen_qform_t* form);
static inline void gen_qform_set(gen_qform_group_t* group, gen_qform_t* R, const gen_qform_t* A);
static inline int gen_qform_equal(gen_qform_group_t* group, const gen_qform_t* A, const gen_qform_t* B);
static inline void gen_qform_inverse(gen_qform_group_t* group, gen_qform_t* form);
static inline void gen_qform_compose(gen_qform_group_t* group, gen_qform_t* R, const gen_qform_t* A, const gen_qform_t* B);
static inline void gen_qform_square(gen_qform_group_t* group, gen_qform_t* R, const gen_qform_t* A);
static inline void gen_qform_cube(gen_qform_group_t* group, gen_qform_t* R, const gen_qform_t* A);
static inline void gen_qform_print(gen_qform_group_t* group, const gen_qform_t* form);

static inline void gen_qform_group_set_discriminant(gen_qform_group_t* this, const mpz_t D);
static inline void gen_qform_reduce(gen_qform_group_t* group, gen_qform_t* form);
static inline int gen_qform_is_ambiguous(gen_qform_group_t* group, const gen_qform_t* form);
static inline int gen_qform_split_ambiguous(gen_qform_group_t* group, mpz_t d, const mpz_t N, const gen_qform_t* form);
static inline int gen_qform_is_primeform(gen_qform_group_t* group, gen_qform_t* form, const int p);

void gen_qform_pow_u32(group_pow_t* pow,
                       gen_qform_t* R,
                       const gen_qform_t* A,
                       uint32_t exp);

/**
 * Inline methods
 */

static inline
void gen_qform_group_init(gen_qform_group_t* this) {
  // initialize group description
  this->desc.group.elem_init = (group_elem_init_f*)&gen_qform_init;
  this->desc.group.elem_clear = (group_elem_clear_f*)&gen_qform_clear;
  this->desc.group.elem_size = sizeof(gen_qform_t);
  this->desc.group.hash32 = (group_hash32_f*)&gen_qform_hash32;
  this->desc.group.set_id = (group_set_id_f*)&gen_qform_set_id;
  this->desc.group.is_id = (group_is_id_f*)&gen_qform_is_id;
  this->desc.group.set = (group_set_f*)&gen_qform_set;
  this->desc.group.equal = (group_equal_f*)&gen_qform_equal;
  this->desc.group.inverse = (group_inverse_f*)&gen_qform_inverse;
  this->desc.group.compose = (group_compose_f*)&gen_qform_compose;
  this->desc.group.square = (group_square_f*)&gen_qform_square;
  this->desc.group.cube = (group_cube_f*)&gen_qform_cube;
  this->desc.group.print = (group_print_f*)&gen_qform_print;
  
  this->desc.discriminant_max_bits = INT_MAX;
  this->desc.set_discriminant = (qform_group_set_discriminant_f*)&gen_qform_group_set_discriminant;
  this->desc.reduce = (qform_reduce_f*)&gen_qform_reduce;
  this->desc.is_primeform = (qform_is_primeform_f*)&gen_qform_is_primeform;
  this->desc.is_ambiguous = (qform_is_ambiguous_f*)&gen_qform_is_ambiguous;
  this->desc.split_ambiguous = (qform_split_ambiguous_f*)&gen_qform_split_ambiguous;
  
  s64_qform_group_init(&this->s64_qform_group);
  s128_qform_group_init(&this->s128_qform_group);
  mpz_qform_group_init(&this->mpz_qform_group);
  
  this->logD = 0;
}

static inline
void gen_qform_group_clear(gen_qform_group_t* this) {
  s64_qform_group_clear(&this->s64_qform_group);
  s128_qform_group_clear(&this->s128_qform_group);
  mpz_qform_group_clear(&this->mpz_qform_group);
  this->logD = 0;
}

static inline
void gen_qform_init(gen_qform_group_t* group, gen_qform_t* this) {
  s64_qform_init(&group->s64_qform_group, &this->s64_qform);
  s128_qform_init(&group->s128_qform_group, &this->s128_qform);
  mpz_qform_init(&group->mpz_qform_group, &this->mpz_qform);
}

static inline
void gen_qform_clear(gen_qform_group_t* group, gen_qform_t* this) {
  s64_qform_clear(&group->s64_qform_group, &this->s64_qform);
  s128_qform_clear(&group->s128_qform_group, &this->s128_qform);
  mpz_qform_clear(&group->mpz_qform_group, &this->mpz_qform);
}

static inline
uint32_t gen_qform_hash32(gen_qform_group_t* group, const gen_qform_t* this) {
  if (group->logD <= s64_qform_group_max_bits) {
    return s64_qform_hash32(&group->s64_qform_group, &this->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    return s128_qform_hash32(&group->s128_qform_group, &this->s128_qform);
  } else {
    return mpz_qform_hash32(&group->mpz_qform_group, &this->mpz_qform);
  }
}

static inline
void gen_qform_set_id(gen_qform_group_t* group, gen_qform_t* this) {
  if (group->logD <= s64_qform_group_max_bits) {
    s64_qform_set_id(&group->s64_qform_group, &this->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    s128_qform_set_id(&group->s128_qform_group, &this->s128_qform);
  } else {
    mpz_qform_set_id(&group->mpz_qform_group, &this->mpz_qform);
  }
}

static inline
int gen_qform_is_id(gen_qform_group_t* group, const gen_qform_t* this) {
  if (group->logD <= s64_qform_group_max_bits) {
    return s64_qform_is_id(&group->s64_qform_group, &this->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    return s128_qform_is_id(&group->s128_qform_group, &this->s128_qform);
  } else {
    return mpz_qform_is_id(&group->mpz_qform_group, &this->mpz_qform);
  }
}

static inline
void gen_qform_set(gen_qform_group_t* group, gen_qform_t* R, const gen_qform_t* A) {
  if (group->logD <= s64_qform_group_max_bits) {
    s64_qform_set(&group->s64_qform_group, &R->s64_qform, &A->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    s128_qform_set(&group->s128_qform_group, &R->s128_qform, &A->s128_qform);
  } else {
    mpz_qform_set(&group->mpz_qform_group, &R->mpz_qform, &A->mpz_qform);
  }
}

static inline
int gen_qform_equal(gen_qform_group_t* group, const gen_qform_t* A, const gen_qform_t* B) {
  if (group->logD <= s64_qform_group_max_bits) {
    return s64_qform_equal(&group->s64_qform_group, &A->s64_qform, &B->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    return s128_qform_equal(&group->s128_qform_group, &A->s128_qform, &B->s128_qform);
  } else {
    return mpz_qform_equal(&group->mpz_qform_group, &A->mpz_qform, &B->mpz_qform);
  }
}

static inline
void gen_qform_inverse(gen_qform_group_t* group, gen_qform_t* this) {
  if (group->logD <= s64_qform_group_max_bits) {
    s64_qform_inverse(&group->s64_qform_group, &this->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    s128_qform_inverse(&group->s128_qform_group, &this->s128_qform);
  } else {
    mpz_qform_inverse(&group->mpz_qform_group, &this->mpz_qform);
  }
}

static inline
void gen_qform_compose(gen_qform_group_t* group, gen_qform_t* R, const gen_qform_t* A, const gen_qform_t* B) {
  if (group->logD <= s64_qform_group_max_bits) {
    s64_qform_compose(&group->s64_qform_group, &R->s64_qform, &A->s64_qform, &B->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    s128_qform_compose(&group->s128_qform_group, &R->s128_qform, &A->s128_qform, &B->s128_qform);
  } else {
    mpz_qform_compose(&group->mpz_qform_group, &R->mpz_qform, &A->mpz_qform, &B->mpz_qform);
  }
}

static inline
void gen_qform_square(gen_qform_group_t* group, gen_qform_t* R, const gen_qform_t* A) {
  if (group->logD <= s64_qform_group_max_bits) {
    s64_qform_square(&group->s64_qform_group, &R->s64_qform, &A->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    s128_qform_square(&group->s128_qform_group, &R->s128_qform, &A->s128_qform);
  } else {
    mpz_qform_square(&group->mpz_qform_group, &R->mpz_qform, &A->mpz_qform);
  }
}

static inline
void gen_qform_cube(gen_qform_group_t* group, gen_qform_t* R, const gen_qform_t* A) {
  if (group->logD <= s64_qform_group_max_bits) {
    s64_qform_cube(&group->s64_qform_group, &R->s64_qform, &A->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    s128_qform_cube(&group->s128_qform_group, &R->s128_qform, &A->s128_qform);
  } else {
    mpz_qform_cube(&group->mpz_qform_group, &R->mpz_qform, &A->mpz_qform);
  }
}

static inline
void gen_qform_print(gen_qform_group_t* group, const gen_qform_t* this) {
  if (group->logD <= s64_qform_group_max_bits) {
    s64_qform_print(&group->s64_qform_group, &this->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    s128_qform_print(&group->s128_qform_group, &this->s128_qform);
  } else {
    mpz_qform_print(&group->mpz_qform_group, &this->mpz_qform);
  }
}

static inline
void gen_qform_group_set_discriminant(gen_qform_group_t* group, const mpz_t D) {
  group->logD = mpz_sizeinbase(D, 2);
  if (group->logD <= s64_qform_group_max_bits) {
    s64_qform_group_set_discriminant(&group->s64_qform_group, D);
  } else if (group->logD <= s128_qform_group_max_bits) {
    s128_qform_group_set_discriminant(&group->s128_qform_group, D);
  } else {
    mpz_qform_group_set_discriminant(&group->mpz_qform_group, D);
  }
}


static inline
void gen_qform_reduce(gen_qform_group_t* group, gen_qform_t* this) {
  if (group->logD <= s64_qform_group_max_bits) {
    s64_qform_reduce(&group->s64_qform_group, &this->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    s128_qform_reduce(&group->s128_qform_group, &this->s128_qform);
  } else {
    mpz_qform_reduce(&group->mpz_qform_group, &this->mpz_qform);
  }
}

static inline
int gen_qform_is_ambiguous(gen_qform_group_t* group, const gen_qform_t* this) {
  if (group->logD <= s64_qform_group_max_bits) {
    return s64_qform_is_ambiguous(&group->s64_qform_group, &this->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    return s128_qform_is_ambiguous(&group->s128_qform_group, &this->s128_qform);
  } else {
    return mpz_qform_is_ambiguous(&group->mpz_qform_group, &this->mpz_qform);
  }
}

static inline
int gen_qform_split_ambiguous(gen_qform_group_t* group, mpz_t d, const mpz_t N, const gen_qform_t* this) {
  if (group->logD <= s64_qform_group_max_bits) {
    return s64_qform_split_ambiguous(&group->s64_qform_group, d, N, &this->s64_qform);
  } else if (group->logD <= s128_qform_group_max_bits) {
    return s128_qform_split_ambiguous(&group->s128_qform_group, d, N, &this->s128_qform);
  } else {
    return mpz_qform_split_ambiguous(&group->mpz_qform_group, d, N, &this->mpz_qform);
  }
}

static inline
int gen_qform_is_primeform(gen_qform_group_t* group, gen_qform_t* this, const int p) {
  if (group->logD <= s64_qform_group_max_bits) {
    return s64_qform_is_primeform(&group->s64_qform_group, &this->s64_qform, p);
  } else if (group->logD <= s128_qform_group_max_bits) {
    return s128_qform_is_primeform(&group->s128_qform_group, &this->s128_qform, p);
  } else {
    return mpz_qform_is_primeform(&group->mpz_qform_group, &this->mpz_qform, p);
  }
}

#endif

