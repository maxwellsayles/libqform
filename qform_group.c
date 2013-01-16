#include "libqform/qform_group.h"

#include "liboptarith/group_pow.h"
#include "liboptarith/math32.h"
#include "liboptarith/primes.h"


/**
 * Return a random unambiguous small prime.
 */
void qform_random_primeform(qform_group_t* group, qform_t* form) {
  int index = rand() % prime_list_count;
  while (!group->is_primeform(group, form, prime_list[index]) ||
	 group->is_ambiguous(group, form)) {
    index = rand() % prime_list_count;
  }
}

/**
 * Compute an unambiguous prime form for a prime with
 * the smallest index >= to 'prime_index'.
 * @return the smallest index greater than or equal to 'prime_index'.
 *         -1 on error.
 */
int qform_next_primeform(qform_group_t* group,
                         qform_t* form,
                         int prime_index) {
  while (!group->is_primeform(group, form, prime_list[prime_index]) ||
	 group->is_ambiguous(group, form)) {
    // Get the next smallest prime.
    prime_index ++;
    if (prime_index >= prime_list_count) {
      return -1;
    }
  }
  return prime_index;
}

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
void qform_pow_u32(group_pow_t* pow,
		   qform_t* R,
		   const qform_t* A,
		   uint32_t exp) {
  qform_group_t* qgroup = (qform_group_t*)pow->group;
  group_t* group = &qgroup->group;

  // Exponentiate by a 16-bit precomputed 2,3 rep using the
  // highest set 16-bits of the exponent.
  const int msb = msb_u32(exp);
  const int b = msb > 15 ? msb+1-16 : 0;
  const uint16_t exp16 = exp >> b;
  const factored_two_three_term16_t* terms = qgroup->pow_reps[exp16];
  const int term_count = qgroup->pow_rep_sizes[exp16];
  group_pow_factored23(pow, pow->E, A, terms, term_count);

  // Exponentiate the rest using binary exponentation.
  if (b > 0) {
    uint32_t m = 1UL << (b-1);
    while (m) {
      group->square(group, pow->E, pow->E);
      if (exp&m) {
	group->compose(group, pow->E, pow->E, A);
      }
      m >>= 1;
    }
  }

  // Copy to output.
  group->set(group, R, pow->E);
}

