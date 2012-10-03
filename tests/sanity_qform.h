#ifndef SANITY_QFORM__INCLUDED
#define SANITY_QFORM__INCLUDED

#include <gmp.h>
#include <limits.h>
#include <stdint.h>

#include "libqform/gen_qform.h"
#include "libqform/mpz_qform.h"
#include "libqform/s128_qform.h"
#include "libqform/s64_qform.h"

typedef struct {
  s64_qform_t s64_qform;
  s128_qform_t s128_qform;
  mpz_qform_t mpz_qform;
  gen_qform_t gen_qform;
} sanity_qform_t;

typedef struct {
  // must be the first member of struct
  qform_group_t desc;
  
  s64_qform_group_t s64_group;
  s128_qform_group_t s128_group;
  mpz_qform_group_t mpz_group;
  gen_qform_group_t gen_group;
  
  int s64_enabled;
  int s128_enabled;

  mpz_t tmp;
  mpz_t tmp2;
} sanity_qform_group_t;

void sanity_qform_group_init(sanity_qform_group_t* this);
void sanity_qform_group_clear(sanity_qform_group_t* this);

static inline void sanity_qform_init(sanity_qform_group_t* group, sanity_qform_t* this);
static inline void sanity_qform_clear(sanity_qform_group_t* group, sanity_qform_t* this);
uint32_t sanity_qform_hash32(sanity_qform_group_t* group, const sanity_qform_t* form);
void sanity_qform_set_id(sanity_qform_group_t* group, sanity_qform_t* form);
int sanity_qform_is_id(sanity_qform_group_t* group, const sanity_qform_t* form);
void sanity_qform_set(sanity_qform_group_t* group, sanity_qform_t* R, const sanity_qform_t* A);
int sanity_qform_equal(sanity_qform_group_t* group, const sanity_qform_t* A, const sanity_qform_t* B);
void sanity_qform_inverse(sanity_qform_group_t* group, sanity_qform_t* form);
void sanity_qform_compose(sanity_qform_group_t* group, sanity_qform_t* R, const sanity_qform_t* A, const sanity_qform_t* B);
void sanity_qform_square(sanity_qform_group_t* group, sanity_qform_t* R, const sanity_qform_t* A);
void sanity_qform_cube(sanity_qform_group_t* group, sanity_qform_t* R, const sanity_qform_t* A);
static inline void sanity_qform_print(sanity_qform_group_t* group, const sanity_qform_t* form);

void sanity_qform_group_set_discriminant(sanity_qform_group_t* this, const mpz_t D);
void sanity_qform_reduce(sanity_qform_group_t* group, sanity_qform_t* form);
int sanity_qform_is_ambiguous(sanity_qform_group_t* group, const sanity_qform_t* form);
int sanity_qform_split_ambiguous(sanity_qform_group_t* group, mpz_t d, const mpz_t N, const sanity_qform_t* form);
int sanity_qform_is_primeform(sanity_qform_group_t* group, sanity_qform_t* form, const int p);

/**
 * Inline methods
 */

static inline void sanity_qform_init(sanity_qform_group_t* group, sanity_qform_t* form) {
  s64_qform_init(&group->s64_group, &form->s64_qform);
  s128_qform_init(&group->s128_group, &form->s128_qform);
  mpz_qform_init(&group->mpz_group, &form->mpz_qform);
  group_t* gen_group = (group_t*)&group->gen_group;
  gen_group->elem_init(gen_group, &form->gen_qform);
}

static inline void sanity_qform_clear(sanity_qform_group_t* group, sanity_qform_t* form) {
  s64_qform_clear(&group->s64_group, &form->s64_qform);
  s128_qform_clear(&group->s128_group, &form->s128_qform);
  mpz_qform_clear(&group->mpz_group, &form->mpz_qform);
  group_t* gen_group = (group_t*)&group->gen_group;
  gen_group->elem_clear(gen_group, &form->gen_qform);
}

static inline void sanity_qform_print(sanity_qform_group_t* group, const sanity_qform_t* form) {
  mpz_qform_print(&group->mpz_group, &form->mpz_qform);
}

#endif // SANTIY_QFORM__INCLUDED


