/*
 * Copyright (c) 2012-2013 Haitham Hassanieh, Piotr Indyk, Dina Katabi,
 *   Eric Price, Massachusetts Institute of Technology
 * Copyright (c) 2012-2013 Jörn Schumacher, ETH Zurich 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */

#ifndef SFFT_H
#define SFFT_H

//#include <google/sparse_hash_map>
#include <tr1/unordered_map>

#include "fftw3.h"
#include "fft.h"
#include "filters.h"

/** Type Definitions **********************************************************/

typedef
std::tr1::unordered_map < int, complex_t, std::tr1::hash < int > > sfft_output;

enum sfft_version
{
  SFFT_VERSION_1,
  SFFT_VERSION_2,
  SFFT_VERSION_3
};

struct sfft_plan
{
  sfft_version version;
  unsigned int n;
  unsigned int k;
  void *data;
};

struct sfft_v1v2_threadlocal_data
{
  /* storage (v1 and v2) */
  int *permute;
  int *permuteb;
  complex_t *x_samp;
  int *score;
  int *hits;
  int *Comb_Approved;
  int *J;
  complex_t *inner_loop_locate_x_sampt;
  real_t *inner_loop_locate_samples;
  real_t **estimate_values_values;

  real_t *nth_element_storage;
  int *nth_int_element_storage;

  /* storage (v2, everything related to Comb filters */
  complex_t *comb_filt_x_sampt;
  real_t *comb_filt_samples;
  std::pair < int, int >*inner_loop_filter_comb_permuted_approved;

  /* fftw plans */
  fftw_plan fftw_plan_comb_filt;
  fftw_plan fftw_plan_inner_loop_locate_location;
  fftw_plan fftw_plan_inner_loop_locate_estimation;

};

struct sfft_v1v2_data
{
  unsigned threads;
  int fftw_opt;

  /* sfft v1.0 and v2.0 */
  int B_loc;
  int B_thresh;
  int B_est;

  int W_Comb;
  int Comb_loops;

  int loops_loc;
  int loops_thresh;
  int loops_est;

  Filter filter;
  Filter filter_est;

  /* storage */
  size_t x_samp_size;
  int *inner_loop_locate_loop_sizes;
  sfft_v1v2_threadlocal_data *threadlocal_data;
};

typedef
sfft_v1v2_data sfft_v1_data;
typedef
sfft_v1v2_data sfft_v2_data;

struct sfft_v3_threadlocal_data
{
  complex_t *gauss_samples;
  complex_t *gauss_perm_samples;
  complex_t *man_samples;
  int *est_freqs;
  complex_t *est_values;
  complex_t *perm_x;
  real_t *nth_element_storage;
  int *nth_int_element_storage;

  fftw_plan fftw_plan_gauss;
  fftw_plan fftw_plan_gauss_perm;
  fftw_plan fftw_plan_man;
};

struct sfft_v3_data
{
  unsigned threads;

  /* sfft v3.0 */
  int B_g1;
  int Gauss_loops;
  int w_g1;
  int B_g2;
  int Gauss2_loops;
  int w_g2;
  int W_Man;
  int Man_loops;
  complex_t *filtert1;
  complex_t *filterf1;
  complex_t *filtert2;
  complex_t *filterf2;

  /* fftw plans */
  int fftw_opt;
  int *fftw_plan_gauss_sizes;
  int *fftw_plan_gauss_perm_sizes;
  int *fftw_plan_man_sizes;

  /* storage */
  sfft_v3_threadlocal_data *threadlocal_data;
};

/** Public API ****************************************************************/

void *sfft_malloc(size_t s);
void sfft_free(void *p);

sfft_plan *sfft_make_plan(int n, int k, sfft_version version,
                          int fftw_optimization);
void sfft_free_plan(sfft_plan * plan);
void sfft_exec(sfft_plan * plan, complex_t * in, sfft_output * out);
void
sfft_exec_many(sfft_plan * plan, int num, complex_t ** in, sfft_output * out);
#endif
