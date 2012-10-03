#include "liboptarith/primes.h"
#include "libqform/qform_group.h"

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
