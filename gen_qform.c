#include <assert.h>

#include "liboptarith/closest_23.h"
#include "liboptarith/group_pow.h"
#include "liboptarith/math32.h"
#include "libqform/gen_qform.h"

/**
 * Exponentiates the high 16 bits using a precomputed 2,3 rep and then
 * exponentiates the remainder using left-to-right binary exponentiation.
 *
 * We use binary rather than a left-to-right NAF since the NAF requires
 * some precomputation when used left-to-right.
 *
 * @param R is the result.
 * @param A is the base.
 * @param exp is the exponent.
 */
void gen_qform_pow_u32(group_pow_t* pow,
                       gen_qform_t* R,
                       const gen_qform_t* A,
                       uint32_t exp) {
  gen_qform_group_t* group = (gen_qform_group_t*)pow->group; 
  const factored_two_three_term16_t* terms = 0;
  int term_count = 0;

  // Exponentiate by a 16-bit precomputed 2,3 rep using the
  // highest set 16-bits of the exponent.
  int msb = msb_u32(exp);
  int b = msb > 15 ? msb+1-16 : 0;
  uint16_t exp16 = exp >> b;
  if (group->logD <= s64_qform_group_max_bits) {
    terms = s64_pow_reps[exp16];
    term_count = s64_pow_rep_sizes[exp16];
  } else if (group->logD <= s128_qform_group_max_bits) {
    terms = s128_pow_reps[exp16];
    term_count = s128_pow_rep_sizes[exp16];
  } else {
    terms = mpz_pow_reps[exp16];
    term_count = mpz_pow_rep_sizes[exp16];
  }
  group_pow_factored23(pow, pow->E, A, terms, term_count);

  // Exponentiate the rest using binary exponentation.
  if (b > 0) {
    uint32_t m = 1UL << (b-1);
    while (m) {
      gen_qform_square(group, pow->E, pow->E);
      if (exp&m) gen_qform_compose(group, pow->E, pow->E, A);
      m >>= 1;
    }
  }

#ifdef _DEBUG
  // Exponentiate using NAF and compare for equality.
  gen_qform_t T;
  gen_qform_init(group, &T);
  group_pow_naf_r2l_u32(pow, &T, A, exp);
  assert(gen_qform_equal(group, &T, pow->E));
  gen_qform_clear(group, &T);
#endif

  // Copy to output.
  gen_qform_set(group, R, pow->E);
}

