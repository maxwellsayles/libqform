#include <stdarg.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

// This comes after so that va_list is defined.
#include <gmp.h>

#include "liboptarith/group.h"
#include "liboptarith/math_mpz.h"
#include "libqform/mpz_qform.h"
#include "libqform/qform_group.h"
#include "libqform/s128_qform.h"
#include "libqform/s64_qform.h"

int rand_seed = 0;

#define verbose 1
#define very_verbose 0

#define cubing_groups 100
#define cubing_ops 100
#define cubing_reps 100

#define min_bits 16
#define max_bits 140
#define total_bits (max_bits+1-min_bits)

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

// gives the time from system on in nanoseconds
static inline uint64_t current_nanos(void) {
#ifdef __linux__
  struct timespec res;
  clock_gettime(CLOCK_MONOTONIC, &res);
  return (res.tv_sec * 1000000000ULL) + res.tv_nsec;
#else
  struct timeval tv;
  gettimeofday(&tv, 0);
  return ((uint64_t)tv.tv_sec * 1000000ULL + (uint64_t)tv.tv_usec) * 1000;
#endif
}

/**
 * Write a gnuplot data file
 */
void write_gnuplot_datfile_u64(const char* filename, const uint64_t* x, const uint64_t* y, int sample_count) {
  int i;
  FILE* f;
  
  f = fopen(filename, "w");
  for (i = 0;  i < sample_count;  i ++) {
    fprintf(f, "%"PRIu64", %"PRIu64"\n", x[i], y[i]);
  }
  fclose(f);
}

typedef struct {
  uint64_t compose_square;
  uint64_t cube;
} cubing_timings_t;

static inline void zero_cubing_timings(cubing_timings_t* this) {
  this->compose_square = 0;
  this->cube = 0;
}

static inline void add_cubing_timings(cubing_timings_t* dst, const cubing_timings_t* src) {
  dst->compose_square += src->compose_square;
  dst->cube += src->cube;
}

void time_cubing_set(cubing_timings_t* timings, int nbits, qform_group_t* qform) {
  group_t* group = &qform->group;
  gmp_randstate_t rands;
  void* form;
  void* tform;
  int i;
  int j;
  mpz_t D;
  uint64_t start;
  uint64_t time;
  
  zero_cubing_timings(timings);
  
  gmp_randinit_default(rands);
  mpz_init(D);
  
  form = group_elem_alloc(group);
  tform = group_elem_alloc(group);
  
  // time cubing
  cprintf("Time for cubing: ");
  fflush(stdout);
  gmp_randseed_ui(rands, rand_seed);
  srand(rand_seed);
  time = 0;
  for (i = 0;  i < cubing_groups;  i ++) {
    if (very_verbose) printf("cube i=%d\n", i);
    
    // generate random semiprime discriminant
    mpz_random_semiprime_discriminant(D, rands, nbits);
    qform->set_discriminant(qform, D);
    
    qform_random_primeform(qform, form);
    start = current_nanos();
    for (j = 0;  j < cubing_ops;  j ++) {
      group->cube(group, form, form);
    }
    time += current_nanos() - start;
  }
  timings->cube = time;
  cprintf("%"PRIu64" us\n", timings->cube);
  
  // time mpz cubing
  cprintf("Time for square+compose: ");
  fflush(stdout);
  gmp_randseed_ui(rands, rand_seed);
  srand(rand_seed);
  time = 0;
  for (i = 0;  i < cubing_groups;  i ++) {
    if (very_verbose) printf("square+compose i=%d\n", i);
    
    // generate random semiprime discriminant
    mpz_random_semiprime_discriminant(D, rands, nbits);
    qform->set_discriminant(qform, D);
    
    qform_random_primeform(qform, form);
    start = current_nanos();
    for (j = 0;  j < cubing_ops;  j ++) {
      group->square(group, tform, form);
      group->compose(group, form, tform, form);
    }
    time += current_nanos() - start;
  }
  timings->compose_square = time;
  cprintf("%"PRIu64" us\n", timings->compose_square);
  
  group_elem_free(group, form);
  group_elem_free(group, tform);
  mpz_clear(D);
}

void time_cubing(void) {
  s64_qform_group_t s64_qform;
  s128_qform_group_t s128_qform;
  mpz_qform_group_t mpz_qform;
  cubing_timings_t timing;
  cubing_timings_t timings[max_bits+1];
  int i;
  int rep;
  uint64_t x[max_bits+1];
  uint64_t y[max_bits+1];
  
  for (i = 0;  i <= max_bits;  i ++) {
    zero_cubing_timings(&timings[i]);
    x[i] = i;
  }
  
  s64_qform_group_init(&s64_qform);
  s128_qform_group_init(&s128_qform);
  mpz_qform_group_init(&mpz_qform);
  
  // load the cpu for a second
  // this primes the os to give us more time slices
  cprintf("Priming the CPU for a while...\n");
  for (i = min_bits;  i <= s64_qform.desc.discriminant_max_bits;  i ++) {
    time_cubing_set(&timing, i, (qform_group_t*)&s64_qform);
  }
  cprintf("\n");
  
  for (rep = 0;  rep < cubing_reps;  rep ++) {
    rand_seed = current_nanos();
    for (i = min_bits;  i <= s64_qform.desc.discriminant_max_bits;  i ++) {
      cprintf("Rep %d on s64_qform for %d bit discriminant.\n", rep, i);
      time_cubing_set(&timing, i, (qform_group_t*)&s64_qform);
      add_cubing_timings(&timings[i], &timing);
      cprintf("\n");
    }
    for (;  i <= s128_qform.desc.discriminant_max_bits;  i ++) {
      cprintf("Rep %d on s128_qform for %d bit discriminant.\n", rep, i);
      time_cubing_set(&timing, i, (qform_group_t*)&s128_qform);
      add_cubing_timings(&timings[i], &timing);
      cprintf("\n");
    }
    for (;  i <= max_bits;  i ++) {
      cprintf("Rep %d on mpz_qform for %d bit discriminant.\n", rep, i);
      time_cubing_set(&timing, i, (qform_group_t*)&mpz_qform);
      add_cubing_timings(&timings[i], &timing);
      cprintf("\n");
    }
  }
  
  // write out cubing times
  for (i = min_bits;  i <= max_bits;  i ++) {
    y[i] = timings[i].cube;
  }
  write_gnuplot_datfile_u64("cube.dat", &x[min_bits], &y[min_bits], total_bits);
  
  // write out compose+square times
  for (i = min_bits;  i <= max_bits;  i ++) {
    y[i] = timings[i].compose_square;
  }
  write_gnuplot_datfile_u64("compose_square.dat", &x[min_bits], &y[min_bits], total_bits);
  
  s64_qform_group_clear(&s64_qform);
  s128_qform_group_clear(&s128_qform);
  mpz_qform_group_clear(&mpz_qform);
}

int main(int argc, char** argv) {
  time_cubing();
  return 0;
}

