#include <gmp.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "libqform/qform_group.h"
#include "libqform/tests/sanity_qform.h"
#include "liboptarith/group_pow.h"
#include "liboptarith/math_mpz.h"

#define iterations 10000
#define groups 100

#define min_bits 16
#define max_bits 140
#define total_bits (max_bits+1-min_bits)

#define min_exp_bits 256
#define max_exp_bits 1024
#define skip_exp_bits 64

#define very_verbose 0

/**
 * Runs qforms from 16bit discriminants up to 120bit discriminants.
 * Randomly chooses compose,square,cube for a random element.
 * If compose is chosen for a destination element, the destination
 * element will be the composition of the other two elements.
 *
 * The test is run with multiple discriminants for each bit sized discriminant.
 */
void test_basic(void) {
  gmp_randstate_t rands;
  sanity_qform_group_t sanity_group;
  sanity_qform_t A;
  sanity_qform_t B;
  sanity_qform_t C;
  mpz_t discriminant;
  int action;
  int dest;
  int i;
  int g;
  int rand_seed;
  int bits;
  
  mpz_init(discriminant);
  
  sanity_qform_group_init(&sanity_group);
  sanity_qform_init(&sanity_group, &A);
  sanity_qform_init(&sanity_group, &B);
  sanity_qform_init(&sanity_group, &C);
  
  rand_seed = time(0);
  srand(rand_seed);
  gmp_randinit_default(rands);
  gmp_randseed_ui(rands, rand_seed);
  
  printf("Performing %d iterations on %d groups for each discriminant bit length.\n", iterations, groups);
  
  for (bits = min_bits;  bits <= max_bits;  bits ++) {
    printf("Using a discriminant of %d bits.\n", bits);
    for (g = 0;  g < groups;  g ++) {
      mpz_random_semiprime_discriminant(discriminant, rands, bits);
      gmp_printf("Discriminant is %Zd\n", discriminant);
      sanity_qform_group_set_discriminant(&sanity_group, discriminant);
      
      // generate three random prime forms for the base
      qform_random_primeform(&sanity_group.desc, &A);
      qform_random_primeform(&sanity_group.desc, &B);
      qform_random_primeform(&sanity_group.desc, &C);
      if (very_verbose) {
	printf("A = ");
	sanity_qform_print(&sanity_group, &A);
	printf("\n");
	printf("B = ");
	sanity_qform_print(&sanity_group, &B);
	printf("\n");
	printf("C = ");
	sanity_qform_print(&sanity_group, &C);
	printf("\n");
      }
      
      // perform random operations of compose, square, and cube on the above
      for (i = 0;  i < iterations;  i ++) {
	action = rand() % 3;
	dest = rand() % 3;
	if (action == 0) {
	  if (dest == 0) {
	    if (very_verbose) printf("A = B*C\n");
	    sanity_qform_compose(&sanity_group, &A, &B, &C);
	  } else if (dest == 1) {
	    if (very_verbose) printf("B = C*A\n");
	    sanity_qform_compose(&sanity_group, &B, &C, &A);
	  } else {
	    if (very_verbose) printf("C = A*B\n");
	    sanity_qform_compose(&sanity_group, &C, &A, &B);
	  }
	} else if (action == 1) {
	  if (dest == 0) {
	    if (very_verbose) printf("A = A^2\n");
	    sanity_qform_square(&sanity_group, &A, &A);
	  } else if (dest == 1) {
	    if (very_verbose) printf("B = B^2\n");
	    sanity_qform_square(&sanity_group, &B, &B);
	  } else {
	    if (very_verbose) printf("C = C^2\n");
	    sanity_qform_square(&sanity_group, &C, &C);
	  }
	} else {
	  if (dest == 0) {
	    if (very_verbose) printf("A = A^3\n");
	    sanity_qform_cube(&sanity_group, &A, &A);
	  } else if (dest == 1) {
	    if (very_verbose) printf("B = B^3\n");
	    sanity_qform_cube(&sanity_group, &B, &B);
	  } else {
	    if (very_verbose) printf("C = C^3\n");
	    sanity_qform_cube(&sanity_group, &C, &C);
	  }
	}
      }
    }
    
    printf("\n");
  }
  
  gmp_randclear(rands);
  sanity_qform_clear(&sanity_group, &A);
  sanity_qform_clear(&sanity_group, &B);
  sanity_qform_clear(&sanity_group, &C);
  sanity_qform_group_clear(&sanity_group);
  
  mpz_clear(discriminant);
  
  printf("All tests completed successfully\n");
}

void test_pow(void) {
  gmp_randstate_t rands;
  sanity_qform_group_t sanity_group;
  sanity_qform_t B;
  sanity_qform_t P1;
  sanity_qform_t P2;
  group_pow_t pow;
  mpz_t discriminant;
  mpz_t ex;
  int exp_bits;
  int rand_seed;
  int bits;
  two_three_term_t* rep;
  int rep_count;
  factored_two_three_term16_t* frep;
  int frep_count;
  
  mpz_init(discriminant);
  mpz_init(ex);
  
  sanity_qform_group_init(&sanity_group);
  sanity_qform_init(&sanity_group, &B);
  sanity_qform_init(&sanity_group, &P1);
  sanity_qform_init(&sanity_group, &P2);
  
  group_pow_init(&pow, &sanity_group.desc.group);
  
  rand_seed = time(0);
  //rand_seed = 0;
  srand(rand_seed);
  gmp_randinit_default(rands);
  gmp_randseed_ui(rands, rand_seed);
  
  printf("Performing %d iterations on %d groups for each discriminant bit length.\n", iterations, groups);
  
  for (bits = min_bits;  bits <= max_bits;  bits ++) {
    printf("Using a discriminant of %d bits.\n", bits);
    // perform exponentiation
    for (exp_bits = min_exp_bits;
	 exp_bits < max_exp_bits;
	 exp_bits += skip_exp_bits) {
      // Generate a random discriminant.
      mpz_random_semiprime_discriminant(discriminant, rands, bits);
      gmp_printf("Discriminant is %Zd\n", discriminant);
      sanity_qform_group_set_discriminant(&sanity_group, discriminant);
      
      // generate large exponents
      mpz_urandomb(ex, rands, exp_bits);
      gmp_printf("ex=%Zd\n", ex);
      
      // generate a random prime form for the base
      qform_random_primeform(&sanity_group.desc, &B);
      printf("B = ");
      sanity_qform_print(&sanity_group, &B);
      printf("\n");
      
      // perform exponentiation by NAF
      group_pow_naf_r2l(&pow, &P1, &B, ex);
      printf("P1 = ");
      sanity_qform_print(&sanity_group, &P1);
      printf("\n");
      
      // perform exponentiation by 2,3 k-closest
      rep = rep_prune_closest(&rep_count, ex, &s64_qform_costs, 512);
      frep = factored_rep(&frep_count, rep, rep_count);
      group_pow_factored23(&pow, &P2, &B, frep, frep_count);
      free(rep);
      free(frep);
      printf("P2 = ");
      sanity_qform_print(&sanity_group, &P2);
      printf("\n");
      
      // Verify forms are equal.
      if (sanity_qform_equal(&sanity_group, &P1, &P2) == 0) {
	printf("Forms disagree!\n");
	exit(-1);
      }
      printf("\n");
    }
    
    printf("\n");
  }
  
  gmp_randclear(rands);
  group_pow_clear(&pow);
  sanity_qform_clear(&sanity_group, &B);
  sanity_qform_clear(&sanity_group, &P1);
  sanity_qform_clear(&sanity_group, &P2);
  sanity_qform_group_clear(&sanity_group);
  
  mpz_clear(ex);
  mpz_clear(discriminant);
  
  printf("All tests completed successfully\n");
}

int main(int argc, char** argv) {
  if (argc != 2) {
    printf("Usage: %s <opt>\n", argv[0]);
    printf("\n");
    printf("Where <opt> is one of:\n");
    printf("  --basic\t test compose, square, cube for 64, 128, mpz, and generic\n");
    printf("  --pow\t\t test NAF against DBNS\n");
    printf("\n");
    return 0;
  }
  
  if (strcmp(argv[1], "--basic") == 0) {
    test_basic();
  } else if (strcmp(argv[1], "--pow") == 0) {
    test_pow();
  } else {
    printf("Unrecognized command line option\n");
  }
  
  return 0;
}

