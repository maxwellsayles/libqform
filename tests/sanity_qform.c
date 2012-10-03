#include <gmp.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "liboptarith/math_mpz.h"
#include "libqform/tests/sanity_qform.h"

static void sanity_qform_dump(sanity_qform_group_t* group, const sanity_qform_t* form, const char* message) {
  printf("Qform INSANITY: %s\n", message);

  if (group->s64_enabled) {
    printf("s64_qform = ");
    s64_qform_print(&group->s64_group, &form->s64_qform);
    printf("\n");
  }
  
  if (group->s128_enabled) {
    printf("s128_qform = ");
    s128_qform_print(&group->s128_group, &form->s128_qform);
    printf("\n");
  }
  
  printf("mpz_qform = ");
  mpz_qform_print(&group->mpz_group, &form->mpz_qform);
  printf("\n");
  
  printf("gen_qform = ");
  gen_qform_print(&group->gen_group, &form->gen_qform);
  printf("\n");
  
  exit(-1);
}

static void sanity_qform_assert(sanity_qform_group_t* group, const sanity_qform_t* form, const char* message) {
  s128_t t;
  s128_t t2;
  char cbuffer[41];
  
  // if s64 is enabled, the generic should be using that type, so compare directly with it
  if (group->gen_group.logD <= s64_qform_group_max_bits) {
    if (!s64_qform_equal(&group->s64_group, &form->s64_qform, &form->gen_qform.s64_qform)) {
      printf("Generic form should be 64bit but does not match s64_qform.\n");
      sanity_qform_dump(group, form, message);
    }
  } else if (group->gen_group.logD <= s128_qform_group_max_bits) {
    if (!s128_qform_equal(&group->s128_group, &form->s128_qform, &form->gen_qform.s128_qform)) {
      printf("Generic form should be 128bit but does not match s128_qform.\n");
      sanity_qform_dump(group, form, message);
    }
  } else if (!mpz_qform_equal(&group->mpz_group, &form->mpz_qform, &form->gen_qform.mpz_qform)) {
    printf("Generic form should be MPZ but does not match mpz_qform.\n");
    sanity_qform_dump(group, form, message);
  }
  
  // test s64_qform and s128_qform against mpz_qform is s64 and s128 are enabled respectively
  if ((group->s64_enabled && (mpz_cmp_s64(form->mpz_qform.a, form->s64_qform.a) != 0 ||
			      mpz_cmp_s64(form->mpz_qform.b, form->s64_qform.b) != 0 ||
			      mpz_cmp_s64(form->mpz_qform.c, form->s64_qform.c) != 0)) ||
      (group->s128_enabled && (mpz_cmp_s64(form->mpz_qform.a, form->s128_qform.a) != 0 ||
			       mpz_cmp_s64(form->mpz_qform.b, form->s128_qform.b) != 0 ||
			       mpz_cmp_s128(form->mpz_qform.c, &form->s128_qform.c) != 0))) {
    
    printf("Forms MISMATCH\n");
    sanity_qform_dump(group, form, message);
  }
  
  // verify that the discriminant matches the group discriminant for s64_qform
  if (group->s64_enabled) {
    // verify s64 discriminant
    mul_s128_s64_s64(&t, form->s64_qform.b, form->s64_qform.b);
    mul_s128_s64_s64(&t2, form->s64_qform.c, form->s64_qform.a);
    shl_s128(&t2);
    shl_s128(&t2);
    sub_s128_s128(&t, &t2);
    if (!is_equal_s128_s64(&t, group->s64_group.D)) {
      printf("s64_group discriminant mismatch\n");
      printf("Should be %"PRId64" but is %"PRId64" instead.\n", group->s64_group.D, t.v0);
      sanity_qform_dump(group, form, message);
    }
  }
  
  // verify that the discriminant matches the group discriminant for s128_qform
  if (group->s128_enabled) {
    // verify s128 discriminant
    mul_s128_s64_s64(&t, form->s128_qform.b, form->s128_qform.b);
    mul_s128_s128_s64(&t2, &form->s128_qform.c, form->s128_qform.a);
    shl_s128(&t2);
    shl_s128(&t2);
    sub_s128_s128(&t, &t2);
    if (!is_equal_s128_s128(&t, &group->s128_group.D)) {
      printf("s64_group discriminant mismatch\n");
      to_decstr_s128(cbuffer, 40, &group->s128_group.D);
      printf("Should be %s", cbuffer);
      to_decstr_s128(cbuffer, 40, &t);
      printf(" but is %s instead.\n", cbuffer);
      sanity_qform_dump(group, form, message);
    }
  }
  
  // verify mpz discriminant
  mpz_mul(group->tmp, form->mpz_qform.b, form->mpz_qform.b);
  mpz_mul(group->tmp2, form->mpz_qform.a, form->mpz_qform.c);
  mpz_mul_2exp(group->tmp2, group->tmp2, 2);
  mpz_sub(group->tmp, group->tmp, group->tmp2);
  if (mpz_cmp(group->tmp, group->mpz_group.D) != 0) {
    printf("mpz_group discriminant mismatch\n");
    gmp_printf("Should be %Zd but is %Zd instead.\n", group->mpz_group.D, group->tmp);
    sanity_qform_dump(group, form, message);
    }
}

/**
 * Initialize a qform group structure
 */
void sanity_qform_group_init(sanity_qform_group_t* this) {
  // initialize group description
  this->desc.group.elem_init = (group_elem_init_f*)&sanity_qform_init;
  this->desc.group.elem_clear = (group_elem_clear_f*)&sanity_qform_clear;
  this->desc.group.elem_size = sizeof(sanity_qform_t);
  this->desc.group.hash32 = (group_hash32_f*)&sanity_qform_hash32;
  this->desc.group.set_id = (group_set_id_f*)&sanity_qform_set_id;
  this->desc.group.is_id = (group_is_id_f*)&sanity_qform_is_id;
  this->desc.group.set = (group_set_f*)&sanity_qform_set;
  this->desc.group.equal = (group_equal_f*)&sanity_qform_equal;
  this->desc.group.inverse = (group_inverse_f*)&sanity_qform_inverse;
  this->desc.group.compose = (group_compose_f*)&sanity_qform_compose;
  this->desc.group.square = (group_square_f*)&sanity_qform_square;
  this->desc.group.cube = (group_cube_f*)&sanity_qform_cube;
  this->desc.group.print = (group_print_f*)&sanity_qform_print;
  
  this->desc.discriminant_max_bits = INT_MAX;
  this->desc.set_discriminant = (qform_group_set_discriminant_f*)&sanity_qform_group_set_discriminant;
  this->desc.reduce = (qform_reduce_f*)&sanity_qform_reduce;
  this->desc.is_primeform = (qform_is_primeform_f*)&sanity_qform_is_primeform;
  this->desc.is_ambiguous = (qform_is_ambiguous_f*)&sanity_qform_is_ambiguous;
  this->desc.split_ambiguous = (qform_split_ambiguous_f*)&sanity_qform_split_ambiguous;
  
  s64_qform_group_init(&this->s64_group);
  s128_qform_group_init(&this->s128_group);
  mpz_qform_group_init(&this->mpz_group);
  gen_qform_group_init(&this->gen_group);
  
  this->s64_enabled = 0;
  this->s128_enabled = 0;
  
  mpz_init(this->tmp);
  mpz_init(this->tmp2);
}

void sanity_qform_group_clear(sanity_qform_group_t* this) {
  mpz_clear(this->tmp);
  mpz_clear(this->tmp2);
  
  s64_qform_group_clear(&this->s64_group);
  s128_qform_group_clear(&this->s128_group);
  mpz_qform_group_clear(&this->mpz_group);
  gen_qform_group_clear(&this->gen_group);
}

void sanity_qform_group_set_discriminant(sanity_qform_group_t* this, const mpz_t D) {
  int nbits = mpz_sizeinbase(D, 2);
  
  this->s64_enabled = 0;
  this->s128_enabled = 0;
  
  if (nbits <= s64_qform_group_max_bits) {
    s64_qform_group_set_discriminant(&this->s64_group, D);
    this->s64_enabled = 1;
  }
  
  if (nbits <= s128_qform_group_max_bits) {
    s128_qform_group_set_discriminant(&this->s128_group, D);
    this->s128_enabled = 1;
  }
  
  mpz_qform_group_set_discriminant(&this->mpz_group, D);
  gen_qform_group_set_discriminant(&this->gen_group, D);
}

uint32_t sanity_qform_hash32(sanity_qform_group_t* group, const sanity_qform_t* form) {
  uint32_t mpz_res = mpz_qform_hash32(&group->mpz_group, &form->mpz_qform);
  uint32_t gen_res = gen_qform_hash32(&group->gen_group, &form->gen_qform);
  if (gen_res != mpz_res) {
    sanity_qform_dump(group, form, "hash32");
  }

  if (group->s64_enabled) {
    uint32_t s64_res = s64_qform_hash32(&group->s64_group, &form->s64_qform);
    if (s64_res != mpz_res) {
      sanity_qform_dump(group, form, "hash32");
    }
  }
  
  if (group->s128_enabled) {
    uint32_t s128_res = s128_qform_hash32(&group->s128_group, &form->s128_qform);
    if (s128_res != mpz_res) {
      sanity_qform_dump(group, form, "hash32");
    }
  }
  
  return mpz_res;
}

void sanity_qform_set_id(sanity_qform_group_t* group, sanity_qform_t* form) {
  if (group->s64_enabled) {
    s64_qform_set_id(&group->s64_group, &form->s64_qform);
  }
  if (group->s128_enabled) {
    s128_qform_set_id(&group->s128_group, &form->s128_qform);
  }
  mpz_qform_set_id(&group->mpz_group, &form->mpz_qform);
  gen_qform_set_id(&group->gen_group, &form->gen_qform);
  
  sanity_qform_assert(group, form, "set_id");
}

int sanity_qform_is_id(sanity_qform_group_t* group, const sanity_qform_t* form) {
  int mpz_res = mpz_qform_is_id(&group->mpz_group, &form->mpz_qform);
  int gen_res = gen_qform_is_id(&group->gen_group, &form->gen_qform);
  if (gen_res != mpz_res) {
    sanity_qform_dump(group, form, "is_id");
  }
  
  if (group->s64_enabled) {
    int s64_res = s64_qform_is_id(&group->s64_group, &form->s64_qform);
    if (s64_res != mpz_res) {
      sanity_qform_dump(group, form, "is_id");
    }
  }
  
  if (group->s128_enabled) {
    int s128_res = s128_qform_is_id(&group->s128_group, &form->s128_qform);
    if (s128_res != mpz_res) {
      sanity_qform_dump(group, form, "is_id");
    }
  }
  
  return mpz_res;
}

void sanity_qform_reduce(sanity_qform_group_t* group, sanity_qform_t* form) {
  if (group->s64_enabled) {
    s64_qform_reduce(&group->s64_group, &form->s64_qform);
  }
  if (group->s128_enabled) {
    s128_qform_reduce(&group->s128_group, &form->s128_qform);
  }
  mpz_qform_reduce(&group->mpz_group, &form->mpz_qform);
  gen_qform_reduce(&group->gen_group, &form->gen_qform);
  
  sanity_qform_assert(group, form, "reduce");
}

int sanity_qform_is_ambiguous(sanity_qform_group_t* group, const sanity_qform_t* form) {
  int mpz_res = mpz_qform_is_ambiguous(&group->mpz_group, &form->mpz_qform);
  int gen_res = gen_qform_is_ambiguous(&group->gen_group, &form->gen_qform);
  if (gen_res != mpz_res) {
    sanity_qform_dump(group, form, "is_ambiguous");
  }
  
  if (group->s64_enabled) {
    int s64_res = s64_qform_is_ambiguous(&group->s64_group, &form->s64_qform);
    if (s64_res != mpz_res) {
      sanity_qform_dump(group, form, "is_ambiguous");
    }
  }
  
  if (group->s128_enabled) {
    int s128_res = s128_qform_is_ambiguous(&group->s128_group, &form->s128_qform);
    if (s128_res != mpz_res) {
      sanity_qform_dump(group, form, "is_ambiguous");
    }
  }
  
  return mpz_res;
}

/**
 * if N is split by the form (a,b), then d is a factor and 1 is returned
 * otherwise 0 is returned.
 */
int sanity_qform_split_ambiguous(sanity_qform_group_t* group, mpz_t d, const mpz_t N, const sanity_qform_t* form) {
  int mpz_res = mpz_qform_split_ambiguous(&group->mpz_group, d, N, &form->mpz_qform);
  int gen_res = gen_qform_split_ambiguous(&group->gen_group, group->tmp, N, &form->gen_qform);
  if (gen_res != mpz_res) {
    sanity_qform_dump(group, form, "split_ambiguous");
  }
  
  if (group->s64_enabled) {
    int s64_res = s64_qform_split_ambiguous(&group->s64_group, group->tmp, N, &form->s64_qform);
    if (mpz_cmp(d, group->tmp) != 0 || s64_res != mpz_res) {
      sanity_qform_dump(group, form, "split_ambiguous");
    }
  }
  
  if (group->s128_enabled) {
    int s128_res = s128_qform_split_ambiguous(&group->s128_group, group->tmp, N, &form->s128_qform);
    if (mpz_cmp(d, group->tmp) != 0 || s128_res != mpz_res) {
      sanity_qform_dump(group, form, "split_ambiguous");
    }
  }
  
  return mpz_res;
}

/**
 * Compute the ideal form for this prime.
 * This is only efficient for small primes.
 */
int sanity_qform_is_primeform(sanity_qform_group_t* group, sanity_qform_t* form, const int p) {
  int mpz_res = mpz_qform_is_primeform(&group->mpz_group, &form->mpz_qform, p);
  int gen_res = gen_qform_is_primeform(&group->gen_group, &form->gen_qform, p);
  if (gen_res != mpz_res) {
    sanity_qform_dump(group, form, "primeform");
  }
  
  if (group->s64_enabled) {
    int s64_res = s64_qform_is_primeform(&group->s64_group, &form->s64_qform, p);
    if (s64_res != mpz_res) {
      sanity_qform_dump(group, form, "primeform");
    }
  }
  
  if (group->s128_enabled) {
    int s128_res = s128_qform_is_primeform(&group->s128_group, &form->s128_qform, p);
    if (s128_res != mpz_res) {
      sanity_qform_dump(group, form, "primeform");
    }
  }
  
  if (mpz_res) {
    sanity_qform_assert(group, form, "primeform");
  }
  
  return mpz_res;
}

void sanity_qform_set(sanity_qform_group_t* group, sanity_qform_t* R, const sanity_qform_t* A) {
  if (group->s64_enabled) {
    s64_qform_set(&group->s64_group, &R->s64_qform, &A->s64_qform);
  }
  
  if (group->s128_enabled) {
    s128_qform_set(&group->s128_group, &R->s128_qform, &A->s128_qform);
  }
  
  mpz_qform_set(&group->mpz_group, &R->mpz_qform, &A->mpz_qform);
  gen_qform_set(&group->gen_group, &R->gen_qform, &A->gen_qform);
  
  sanity_qform_assert(group, R, "set");
}

int sanity_qform_equal(sanity_qform_group_t* group, const sanity_qform_t* A, const sanity_qform_t* B) {
  int mpz_res = mpz_qform_equal(&group->mpz_group, &A->mpz_qform, &B->mpz_qform);
  int gen_res = gen_qform_equal(&group->gen_group, &A->gen_qform, &B->gen_qform);
  if (gen_res != mpz_res) {
    printf("MISMATCH on gen_qform_equal\n");
    exit(-1);
  }
  
  if (group->s64_enabled) {
    int s64_res = s64_qform_equal(&group->s64_group, &A->s64_qform, &B->s64_qform);
    if (s64_res != mpz_res) {
      printf("MISMATCH on s64_qform_equal\n");
      exit(-1);
    }
  }
  
  if (group->s128_enabled) {
    int s128_res = s128_qform_equal(&group->s128_group, &A->s128_qform, &B->s128_qform);
    if (s128_res != mpz_res) {
      printf("MISMATCH on s128_qform_equal\n");
      exit(-1);
    }
  }
  
  return mpz_res;
}

void sanity_qform_inverse(sanity_qform_group_t* group, sanity_qform_t* form) {   
  if (group->s64_enabled) {
    s64_qform_inverse(&group->s64_group, &form->s64_qform);
  }
  if (group->s128_enabled) {
    s128_qform_inverse(&group->s128_group, &form->s128_qform);
  }
  mpz_qform_inverse(&group->mpz_group, &form->mpz_qform);
  gen_qform_inverse(&group->gen_group, &form->gen_qform);
  
  sanity_qform_assert(group, form, "inverse");
}

/**
 * Computes a reduced ideal equivalent to the product of two ideals
 * using the NUCOMP algorithm of Shanks
 */
void sanity_qform_compose(sanity_qform_group_t* group, sanity_qform_t* R, const sanity_qform_t* A, const sanity_qform_t* B) {
  if (group->s64_enabled) {
    s64_qform_compose(&group->s64_group, &R->s64_qform, &A->s64_qform, &B->s64_qform);
  }
  if (group->s128_enabled) {
    s128_qform_compose(&group->s128_group, &R->s128_qform, &A->s128_qform, &B->s128_qform);
  }
  mpz_qform_compose(&group->mpz_group, &R->mpz_qform, &A->mpz_qform, &B->mpz_qform);
  gen_qform_compose(&group->gen_group, &R->gen_qform, &A->gen_qform, &B->gen_qform);
  
  sanity_qform_assert(group, R, "compose");
}

/**
 * Computes a reduced ideal equivalent to the square of an ideal
 * using the NUCOMP algorithm of Shanks.
 */
void sanity_qform_square(sanity_qform_group_t* group, sanity_qform_t* R, const sanity_qform_t* A) {
  if (group->s64_enabled) {
    s64_qform_square(&group->s64_group, &R->s64_qform, &A->s64_qform);
  }
  if (group->s128_enabled) {
    s128_qform_square(&group->s128_group, &R->s128_qform, &A->s128_qform);
  }
  mpz_qform_square(&group->mpz_group, &R->mpz_qform, &A->mpz_qform);
  gen_qform_square(&group->gen_group, &R->gen_qform, &A->gen_qform);
  
  sanity_qform_assert(group, R, "square");
}

/**
 * Computes a reduced ideal equivalent to the cube of an ideal
 * using an adaptation of Shanks' NUCOMP algorithm.
 */
void sanity_qform_cube(sanity_qform_group_t* group, sanity_qform_t* R, const sanity_qform_t* A) {
  if (group->s64_enabled) {
    s64_qform_cube(&group->s64_group, &R->s64_qform, &A->s64_qform);
  }
  if (group->s128_enabled) {
    s128_qform_cube(&group->s128_group, &R->s128_qform, &A->s128_qform);
  }
  mpz_qform_cube(&group->mpz_group, &R->mpz_qform, &A->mpz_qform);
  gen_qform_cube(&group->gen_group, &R->gen_qform, &A->gen_qform);
  
  sanity_qform_assert(group, R, "cube");
}


