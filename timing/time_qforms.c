#if !defined(NO_PARI)
#include <pari/pari.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

// This comes here so that va_list is defined.
#include <gmp.h>

#include "liboptarith/group.h"
#include "liboptarith/math_mpz.h"
#include "libqform/mpz_qform.h"
#include "libqform/qform_group.h"
#include "libqform/s128_qform.h"
#include "libqform/s64_qform.h"

#define verbose 1
#define very_verbose 0

#define time_s64 0
#define time_s128 1
#define time_mpz 0
#define time_pari 0

#define qform_groups 1000
#define qform_ops 1000
#define qform_reps 1

#define min_bits 16
#define max_bits 140
#define total_bits (max_bits+1-min_bits)

#define min(a,b) (((a)<(b)) ? (a) : (b))

int cprintf(const char* fmt, ...) {
  int res = 0;
  va_list args;
  if (verbose) {
    va_start(args, fmt);
    res = gmp_vprintf(fmt, args);
    va_end(args);
  }
  return res;
}

/// Convert the mpz_t to a string and then to a GEN.
/// This is inefficient, but works, and isn't timed anyways.
#if !defined(NO_PARI)
GEN to_gen(mpz_t x) {
#if GMP_LIMB_BITS == 64
  if (x->_mp_size < 0) {
    printf("Does not support negative integers.\n");
    exit(-1);
  } else if (x->_mp_size == 0) {
    return stoi(0);
  } else if (x->_mp_size == 1) {
    return stoi(x->_mp_d[0]);
  } else if (x->_mp_size == 2) {
    return mkintn(4,
		  (uint64_t)(x->_mp_d[1] >> 32) & 0xFFFFFFFF,
		  (uint64_t)x->_mp_d[1] & 0xFFFFFFFF,
		  (uint64_t)(x->_mp_d[0] >> 32) & 0xFFFFFFFF,
		  (uint64_t)x->_mp_d[0] & 0xFFFFFFFF);
  } else if (x->_mp_size == 3) {
    return mkintn(6,
		  (uint64_t)(x->_mp_d[2] >> 32) & 0xFFFFFFFF,
		  (uint64_t)x->_mp_d[2] & 0xFFFFFFFF,
		  (uint64_t)(x->_mp_d[1] >> 32) & 0xFFFFFFFF,
		  (uint64_t)x->_mp_d[1] & 0xFFFFFFFF,
		  (uint64_t)(x->_mp_d[0] >> 32) & 0xFFFFFFFF,
		  (uint64_t)x->_mp_d[0] & 0xFFFFFFFF);
  } else {
    printf("Unsupported integer size.\n");
    exit(-1);
  }
#else
  // Hack for non 64-bit Intel platforms.
  char c[1024];
  mpz_get_str(c, 10, x);
  GEN y = strtoi(c);
  return y;
#endif
}
#endif

// gives the time from system on in nanoseconds
static inline uint64_t current_nanos(void) {
#ifdef __linux__
  struct timespec res;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &res);
  return (res.tv_sec * 1000000000ULL) + res.tv_nsec;
#else
  struct timeval tv;
  gettimeofday(&tv, 0);
  return ((uint64_t)tv.tv_sec * 1000000ULL + (uint64_t)tv.tv_usec) * 1000;
#endif
}

/// Write a gnuplot data file
void write_gnuplot_datfile_u64(const char* filename,
			       const uint64_t* x,
			       const uint64_t* y,
			       int sample_count) {
  int i;
  FILE* f;
  f = fopen(filename, "w");
  for (i = 0;  i < sample_count;  i ++) {
    fprintf(f, "%"PRIu64", %.5f\n",
	    x[i],
	    (double)y[i] / (double)(qform_groups * qform_ops * qform_reps));
  }
  fclose(f);
}

typedef struct {
  uint64_t compose;
  uint64_t square;
  uint64_t cube;
} qform_timings_t;

static inline void zero_qform_timings(qform_timings_t* this) {
  this->compose = 0;
  this->square = 0;
  this->cube = 0;
}

void time_qform_set(qform_timings_t* timings,
		    qform_group_t* qform_group,
		    const int nbits,
		    const uint64_t rand_seed,
		    const char* name) {
  group_t* group = &qform_group->group;
  gmp_randstate_t rands;
  int i;
  int j;
  mpz_t D;
  uint64_t time;
  uint64_t start;
  void* A;
  void* B;
  void* C;

  gmp_randinit_default(rands);
  mpz_init(D);

  A = group_elem_alloc(group);
  B = group_elem_alloc(group);
  C = group_elem_alloc(group);

  cprintf("Running qform timings for %d bit discriminants.\n", nbits);

  // time compose
  cprintf("Time for %d-bit compose-%s: ", nbits, name);
  fflush(stdout);
  gmp_randseed_ui(rands, rand_seed);
  srand(rand_seed);
  time = 0;
  for (i = 0;  i < qform_groups;  i ++) {
    if (very_verbose) cprintf("%dbit compose i=%d\n", nbits, i);

    // generate random semiprime discriminant
    mpz_random_semiprime_discriminant(D, rands, nbits);
    qform_group->set_discriminant(qform_group, D);

    qform_random_primeform(qform_group, A);
    group->set(group, B, A);
    start = current_nanos();
    for (j = 0;  j < qform_ops;  j ++) {
      group->compose(group, C, B, A);
      group->set(group, A, B);
      group->set(group, B, C);
    }
    time += current_nanos() - start;
  }
  timings->compose += time;
  cprintf("%"PRIu64" us\n", time);

  // time square
  cprintf("Time for %d-bit square-%s: ", nbits, name);
  fflush(stdout);
  gmp_randseed_ui(rands, rand_seed);
  srand(rand_seed);
  time = 0;
  for (i = 0;  i < qform_groups;  i ++) {
    if (very_verbose) cprintf("%dbit square i=%d\n", nbits, i);

    // generate random semiprime discriminant
    mpz_random_semiprime_discriminant(D, rands, nbits);
    qform_group->set_discriminant(qform_group, D);

    qform_random_primeform(qform_group, A);
    start = current_nanos();
    for (j = 0;  j < qform_ops;  j ++) {
      group->square(group, A, A);
    }
    time += current_nanos() - start;
  }
  timings->square += time;
  cprintf("%"PRIu64" us\n", time);

  // time cube
  cprintf("Time for %d-bit cube-%s: ", nbits, name);
  fflush(stdout);
  gmp_randseed_ui(rands, rand_seed);
  srand(rand_seed);
  time = 0;
  for (i = 0;  i < qform_groups;  i ++) {
    if (very_verbose) cprintf("%dbits cube i=%d\n", nbits, i);

    // generate random semiprime discriminant
    mpz_random_semiprime_discriminant(D, rands, nbits);
    qform_group->set_discriminant(qform_group, D);

    qform_random_primeform(qform_group, A);
    start = current_nanos();
    for (j = 0;  j < qform_ops;  j ++) {
      group->cube(group, A, A);
    }
    time += current_nanos() - start;
  }
  timings->cube += time;
  cprintf("%"PRIu64" us\n", time);

  group_elem_free(group, A);
  group_elem_free(group, B);
  group_elem_free(group, C);

  mpz_clear(D);
}

#if !defined(NO_PARI)
void time_pari_set(qform_timings_t* timings,
		   const int nbits,
		   const uint64_t rand_seed,
		   const char* name) {
  gmp_randstate_t rands;
  int i;
  int j;
  mpz_t D;
  uint64_t time;
  uint64_t start;
  GEN A;
  GEN B;
  GEN C;

  mpz_qform_group_t qform_group;
  mpz_qform_t P;

  mpz_qform_group_init(&qform_group);
  mpz_qform_init(&qform_group, &P);

  gmp_randinit_default(rands);
  mpz_init(D);

  cprintf("Running qform timings for %d bit discriminants.\n", nbits);

  // time compose
  cprintf("Time for %dbit compose-%s: ", nbits, name);
  fflush(stdout);
  gmp_randseed_ui(rands, rand_seed);
  srand(rand_seed);
  time = 0;
  for (i = 0;  i < qform_groups;  i ++) {
    if (very_verbose) cprintf("%dbit compose i=%d\n", nbits, i);

    // Generate random semiprime discriminant and primeform.
    mpz_random_semiprime_discriminant(D, rands, nbits);
    mpz_qform_group_set_discriminant(&qform_group, D);
    qform_random_primeform(&qform_group.desc, &P);

    // Create Pari objects.
    pari_sp ltop = avma;
    GEN a = to_gen(P.a);
    GEN b = to_gen(P.b);
    GEN c = to_gen(P.c);
    A = qfi(a, b, c);
    B = qfi(a, b, c);

    start = current_nanos();
    for (j = 0;  j < qform_ops;  j ++) {
      C = gmul(B, A);
      A = B;
      B = C;
    }
    time += current_nanos() - start;
    avma = ltop;  // Clean up Pari's stack.
  }
  timings->compose += time;
  cprintf("%"PRIu64" us\n", time);

  // time square
  cprintf("Time for %dbit square-%s: ", nbits, name);
  fflush(stdout);
  gmp_randseed_ui(rands, rand_seed);
  srand(rand_seed);
  time = 0;
  for (i = 0;  i < qform_groups;  i ++) {
    if (very_verbose) cprintf("%dbit square i=%d\n", nbits, i);

    // Generate random semiprime discriminant and primeform.
    mpz_random_semiprime_discriminant(D, rands, nbits);
    mpz_qform_group_set_discriminant(&qform_group, D);
    qform_random_primeform(&qform_group.desc, &P);

    // Create Pari objects.
    pari_sp ltop = avma;
    GEN a = to_gen(P.a);
    GEN b = to_gen(P.b);
    GEN c = to_gen(P.c);
    A = qfi(a, b, c);
    
    start = current_nanos();
    for (j = 0;  j < qform_ops;  j ++) {
      A = gsqr(A);
    }
    time += current_nanos() - start;
    avma = ltop;
  }
  timings->square += time;
  cprintf("%"PRIu64" us\n", time);

  // time cube
  cprintf("Time for %dbit cube-%s: ", nbits, name);
  fflush(stdout);
  gmp_randseed_ui(rands, rand_seed);
  srand(rand_seed);
  time = 0;
  for (i = 0;  i < qform_groups;  i ++) {
    if (very_verbose) cprintf("%dbits cube i=%d\n", nbits, i);

    // Generate random semiprime discriminant and primeform.
    mpz_random_semiprime_discriminant(D, rands, nbits);
    mpz_qform_group_set_discriminant(&qform_group, D);
    qform_random_primeform(&qform_group.desc, &P);

    // Create Pari objects.
    pari_sp ltop = avma;
    GEN a = to_gen(P.a);
    GEN b = to_gen(P.b);
    GEN c = to_gen(P.c);
    A = qfi(a, b, c);

    start = current_nanos();
    for (j = 0;  j < qform_ops;  j ++) {
      A = gpowgs(A, 3);
    }
    time += current_nanos() - start;
    avma = ltop;
  }
  timings->cube += time;
  cprintf("%"PRIu64" us\n", time);

  mpz_qform_clear(&qform_group, &P);
  mpz_qform_group_clear(&qform_group);
  mpz_clear(D);
}
#endif

void output_timings(const qform_timings_t* timings, const int j,
		    const char* compose_file,
		    const char* square_file,
		    const char* cube_file) {
  int i;
  uint64_t x[max_bits+1];
  uint64_t y[max_bits+1];
  for (i = min_bits; i <= j; i++) {
    x[i] = i;
    y[i] = timings[i].compose;
  }
  write_gnuplot_datfile_u64(compose_file,
			    &x[min_bits], &y[min_bits], j - min_bits + 1);

  for (i = min_bits; i <= j; i++) {
    y[i] = timings[i].square;
  }
  write_gnuplot_datfile_u64(square_file,
			    &x[min_bits], &y[min_bits], j - min_bits + 1);

  for (i = min_bits; i <= j; i++) {
    y[i] = timings[i].cube;
  }
  write_gnuplot_datfile_u64(cube_file,
			    &x[min_bits], &y[min_bits], j - min_bits + 1);
}

void time_qforms(void) {
  qform_timings_t timings_s64[max_bits+1];
  qform_timings_t timings_s128[max_bits+1];
  qform_timings_t timings_mpz[max_bits+1];
  qform_timings_t timings_pari[max_bits+1];
  s64_qform_group_t s64_qform;
  s128_qform_group_t s128_qform;
  mpz_qform_group_t mpz_qform;
  int i;
  int j;
  int rep;
  uint64_t rand_seed;

  s64_qform_group_init(&s64_qform);
  s128_qform_group_init(&s128_qform);
  mpz_qform_group_init(&mpz_qform);

  // Prime the CPU for 3 seconds
  printf("Priming CPU for 3 seconds.\n");
  uint64_t start = current_nanos();
  while (current_nanos() - start < 3000000000) {
  }

  // Zero timings.
  for (i = 0;  i <= max_bits;  i ++) {
    zero_qform_timings(&timings_s64[i]);
    zero_qform_timings(&timings_s128[i]);
    zero_qform_timings(&timings_mpz[i]);
    zero_qform_timings(&timings_pari[i]);
  }

  for (rep = 0;  rep < qform_reps;  rep ++) {
    rand_seed = rep;

#if (time_s64 == 1)
    // Run set using s64 implementation.
    j = min(s64_qform.desc.discriminant_max_bits, max_bits);
    for (i = min_bits; i <= j; i++) {
      cprintf("Rep %d on s64_qform for %d bit discriminant.\n", rep, i);
      time_qform_set(&timings_s64[i], (qform_group_t*)&s64_qform,
		     i, rand_seed, "64");
      cprintf("\n");
    }
#endif

#if (time_s128 == 1)
    // Run set using s128 implementation.
    j = min(s128_qform.desc.discriminant_max_bits, max_bits);
    for (i = min_bits; i <= j; i++) {
      cprintf("Rep %d on s128_qform for %d bit discriminant.\n", rep, i);
      time_qform_set(&timings_s128[i], (qform_group_t*)&s128_qform,
		     i, rand_seed, "128");
      cprintf("\n");
    }
#endif

#if (time_mpz == 1)
    // Run set using MPZ implementation.
    for (i = min_bits; i <= max_bits; i++) {
      cprintf("Rep %d on mpz_qform for %d bit discriminant.\n", rep, i);
      time_qform_set(&timings_mpz[i], (qform_group_t*)&mpz_qform,
		     i, rand_seed, "mpz");
      cprintf("\n");
    }
#endif

#if (time_pari == 1)
    // Run set using Pari implementation.
    for (i = min_bits; i <= max_bits; i++) {
      cprintf("Rep %d on Pari for %d bit discriminant.\n", rep, i);
      time_pari_set(&timings_pari[i], i, rand_seed, "pari");
      cprintf("\n");
    }
#endif
  }

  // output data
#if (time_s64 == 1)
  j = min(s64_qform.desc.discriminant_max_bits, max_bits);
  output_timings(timings_s64, j,
		 "compose-64.dat", "square-64.dat", "cube-64.dat");
#endif
#if (time_s128 == 1)
  j = min(s128_qform.desc.discriminant_max_bits, max_bits);
  output_timings(timings_s128, j,
		 "compose-128.dat", "square-128.dat", "cube-128.dat");
#endif
#if (time_mpz == 1)
  output_timings(timings_mpz, max_bits,
		 "compose-mpz.dat", "square-mpz.dat", "cube-mpz.dat");
#endif
#if (time_pari == 1)
  output_timings(timings_pari, max_bits,
		 "compose-pari.dat", "square-pari.dat", "cube-pari.dat");
#endif

  s64_qform_group_clear(&s64_qform);
  s128_qform_group_clear(&s128_qform);
  mpz_qform_group_clear(&mpz_qform);
}

int main(int argc, char** argv) {
#if (time_pari == 1)
  pari_init(1<<29, 1<<20);
#endif

  time_qforms();
  return 0;
}

