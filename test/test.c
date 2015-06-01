#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <sys/types.h>
#include <unistd.h>
#include "sfft.h"
#include "fftw.h"

#define PI 3.14159265

int main (void)
{
  sfft_plan * myPlan;
  int i;
  int n = 8192;
  int k = 50;
  myPlan = sfft_make_plan (n, k, SFFT_VERSION_2, FFTW_ESTIMATE);

  /* input section */
  complex_t * input_vector = (complex_t *)sfft_malloc(n);
  printf("old: %p\n", input_vector);
  srand(17);
  srand48( time(NULL) ^ (getpid() * 171717));
  int * LARGE_FREQ = (int *) malloc (k * sizeof(*LARGE_FREQ));
  complex_t * x_f = (complex_t *) calloc (n, sizeof(*x_f));
  for (i = 0; i < k; i++) {
    LARGE_FREQ[i] = (unsigned) floor( drand48() * n);
    x_f[LARGE_FREQ[i]] = 1.0;
  }
  fftw_dft(input_vector, n, x_f, 1);
  printf("new: %p\n", input_vector);

  sfft_output output_vector;
  sfft_exec(myPlan, input_vector, &output_vector);
  printf("new2: %p\n", input_vector);

//  sfft_free(input_vector);
  sfft_free_plan(myPlan);

  return 0;
}
