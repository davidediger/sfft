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

#include <cstring>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <cmath>

#include <omp.h>

#include "intrinsics.h"
#include "computefourier-1.0-2.0.h"
#include "filters.h"
#include "utils.h"
#include "timer.h"

#include "flopcount.h"
#include "profiling_tools.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

inline int timesmod(const int &x, const int &a, const int &n)
{
  return int ((((long long int)x) * a) % n);
}

int
Comb_Filt(sfft_v1v2_data * data, complex_t * origx, int n, int num,
          int W_Comb, int *Comb_Approved)
{

  assert(n % W_Comb == 0);

  unsigned tid = omp_get_thread_num();
  sfft_v1v2_threadlocal_data *tldata = &data->threadlocal_data[tid];

  complex_t *x_sampt = tldata->comb_filt_x_sampt;
  real_t *samples = tldata->comb_filt_samples;

  int sigma = n / W_Comb;
  int offset = (unsigned)floor(drand48() * sigma);

  for (int i = 0; i < W_Comb; i++)
    {
      x_sampt[i] = origx[offset + i * sigma];
    }

  fftw_execute(tldata->fftw_plan_comb_filt);
  FLOPCOUNT_INCREMENT_FFTW_PLAN(tldata->fftw_plan_comb_filt);

  for (int i = 0; i < W_Comb; i++)
    {
      samples[i] = cabs2(x_sampt[i]);
      FLOPCOUNT_INCREMENT(2 + 1);
    }

  find_largest_indices(Comb_Approved, num, samples, W_Comb,
                       tldata->nth_element_storage);

  return 0;
}

/*
Find indices that map to J , i.e., lie within n/(2B) of (J * n/B) after
permutation.

For each such i, increment score[i] and append to hits if score[i] reaches
loop_threshold.
*/
int
inner_loop_filter_regular(int *J, int n, int num, int B, int a, int ai,
                          int b, int loop_threshold, int *score, int *hits,
                          int &hits_found)
{
  // Given the set of large samples, find the locations in [n] that map there
  // and output them

  for (int i = 0; i < num; i++)
    {
      int low, high;
      low = (int (ceil((J[i] - 0.5) * n / B)) + n)%n;
      high = (int (ceil((J[i] + 0.5) * n / B)) + n)%n;
      FLOPCOUNT_INCREMENT(2 * 3);
      int loc = timesmod(low, a, n);
      for (int j = low; j != high; j = (j + 1) % n)
        {
          score[loc]++;
          if (score[loc] == loop_threshold)
            hits[hits_found++] = loc;
          loc = (loc + a) % n;
        }
    }

  return 0;
}

/*
Find indices that (1) map to J under the permutation and (2) lie in
Comb_Approved mod W_Comb.

For each such i, increment hits[i] and append to hits_found if hits[i]
reaches loop_threshold.
*/
int
inner_loop_filter_Comb(sfft_v1v2_data * data, int *J, int n, int num,
                       int B, int a, int ai, int b, int loop_threshold,
                       int *score, int *hits, int &hits_found,
                       int *Comb_Approved, int num_Comb, int W_Comb)
{
  unsigned tid = omp_get_thread_num();
  sfft_v1v2_threadlocal_data *tldata = &data->threadlocal_data[tid];

  std::pair < int, int >*permuted_approved =
    tldata->inner_loop_filter_comb_permuted_approved;

  for (int m = 0; m < num_Comb; m++)
    {
      int prev = timesmod(Comb_Approved[m], ai, W_Comb);
      permuted_approved[m] =
        std::make_pair(prev, timesmod(prev, a, n));
    }
  std::sort(permuted_approved, permuted_approved + num_Comb);

  // compute intersection of permuted_approved and indices close to J * n/B, then invert to get true locations.

  for (int i = 0; i < num; i++)
    {
      int low, high;
      low = (int (ceil((J[i] - 0.5) * n / B)) + n)%n;
      high = (int (ceil((J[i] + 0.5) * n / B)) + n)%n;
      FLOPCOUNT_INCREMENT(2 * 3);
      int index =
        int (std::upper_bound
             (permuted_approved, permuted_approved + num_Comb,
              std::make_pair(low % W_Comb,
                             -1)) - permuted_approved);
      int location = low - (low % W_Comb);
      int locinv = timesmod(location, a, n);
      for (int j = index;; j++)
        {
          if (j == num_Comb)
            {
              j -= num_Comb;
              location = (location + W_Comb) % n;
              locinv = timesmod(location, a, n);
            }
          int approved_loc =
            location + permuted_approved[j].first;
          if ((low < high
               && (approved_loc >= high || approved_loc < low))
              || (low > high
                  && (approved_loc >= high
                      && approved_loc < low)))
            break;
          int loc = (locinv + permuted_approved[j].second) % n;
          score[loc]++;
          if (score[loc] == loop_threshold)
            hits[hits_found++] = loc;
        }
    }

  return 0;
}

/*
  Inner loop of the algorithm, part one.

  n-dimensional origx
  permute the fourier spectrum, take the first w coordinates.
  dot with the filter
  B-dimensional FFT
  return the top num samples.
 */
int
inner_loop_locate(sfft_v1v2_data * data, complex_t * origx, int n,
                  const Filter & loc_filter, const Filter & est_filter,
                  int num, int B_loc, int B_est, int *a_vec, int *ai_vec,
                  int *b_vec, complex_t * x_samp, int *J, int loops,
                  int location_loops, int loop_threshold, int *score,
                  int *hits, int &hits_found, int *Comb_Approved,
                  int num_Comb, int W_Comb)
{
  unsigned tid = omp_get_thread_num();
  sfft_v1v2_threadlocal_data *tldata = &data->threadlocal_data[tid];

  complex_t *x_sampt = tldata->inner_loop_locate_x_sampt;
  real_t *samples = tldata->inner_loop_locate_samples;
  size_t x_sampt_size = data->x_samp_size;

  double *d_x_sampt = (double *)x_sampt;
  double *d_orig_x = (double *)origx;
  memset(d_x_sampt, 0, x_sampt_size * sizeof(x_sampt[0]));

  assert(n % B_loc == 0);
  assert(n % B_est == 0);

  const __m128d signs = _mm_set_pd(-1.0, 1.0);
  const int n2_m_1 = 2 * n - 1;	// can be used for fast "mod 2n" calculations

  /* Permutation and filter application */
  for (int j = 0; j < loops; j++)
    {
      int perform_location = (j < location_loops);
      Filter filter = perform_location ? loc_filter : est_filter;
      int B = perform_location ? B_loc : B_est;
      unsigned int B2_m_1 = 2 * B - 1;
      int offset =
        2 * (MIN(j, location_loops) * B_loc +
             MAX(0, j - location_loops) * B_est);

      double *d_filter = (double *)filter.time;

      int ai = ai_vec[j];
      int b = b_vec[j];

      //Permute, dot, collate all in one loop.
      unsigned int index = b;
      int size = 2 * filter.sizet;

      for (int i = 0; i < size;)
        {
          __m128d ab = _mm_load_pd(d_orig_x + index);
          __m128d cd = _mm_load_pd(d_filter + i);
          __m128d dc = _mm_shuffle_pd(cd, cd, 1);

          __m128d ac_bd = _mm_mul_pd(ab, cd);
          __m128d ad_bc = _mm_mul_pd(ab, dc);
          __m128d ac_mbd = _mm_mul_pd(ac_bd, signs);

          __m128d ab_times_cd = _mm_hadd_pd(ac_mbd, ad_bc);

          unsigned int i_mod_B_p_offset = (i & B2_m_1) + offset;
          __m128d xy = _mm_load_pd(d_x_sampt + i_mod_B_p_offset);
          __m128d st = _mm_add_pd(xy, ab_times_cd);
          _mm_store_pd(d_x_sampt + i_mod_B_p_offset, st);

          index = (index + 2 * ai) & n2_m_1;
          i += 2;
        }
    }

  FLOPCOUNT_INCREMENT((location_loops * loc_filter.sizet +
                       (loops -
                        location_loops) * est_filter.sizet) *
                      (4 /*mults */  + 4 /*adds/subs */ ));

  /* FFTs */

  fftw_execute(tldata->fftw_plan_inner_loop_locate_location);
  FLOPCOUNT_INCREMENT_FFTW_PLAN
  (tldata->fftw_plan_inner_loop_locate_location);
  fftw_execute(tldata->fftw_plan_inner_loop_locate_estimation);
  FLOPCOUNT_INCREMENT_FFTW_PLAN
  (tldata->fftw_plan_inner_loop_locate_estimation);

  assert(x_sampt_size % 4 == 0);
  for (unsigned j = 0; j < x_sampt_size; j += 2)
    {
      __m128d ab = _mm_load_pd((double *)(x_samp + j));
      __m128d cd = _mm_load_pd((double *)(x_samp + j + 1));

      __m128d ab_square = _mm_mul_pd(ab, ab);
      __m128d cd_square = _mm_mul_pd(cd, cd);

      __m128d r = _mm_hadd_pd(ab_square, cd_square);

      _mm_store_pd(samples + j, r);
    }
  FLOPCOUNT_INCREMENT(x_sampt_size * (2 + 1));

  for (int j = 0; j < loops; j++)
    {
      int perform_location = (j < location_loops);
      int B = perform_location ? B_loc : B_est;
      int offset = MIN(j, location_loops) * B_loc + MAX(0,
                   j -
                   location_loops)
                   * B_est;

      find_largest_indices(J, num, samples + offset, B,
                           tldata->nth_element_storage);

      if (perform_location)
        {
          if (!WITH_COMB)
            {
              inner_loop_filter_regular(J, n, num, B,
                                        a_vec[j], ai_vec[j],
                                        b_vec[j],
                                        loop_threshold, score,
                                        hits, hits_found);
            }
          else
            {
              inner_loop_filter_Comb(data, J, n, num, B,
                                     a_vec[j], ai_vec[j],
                                     b_vec[j], loop_threshold,
                                     score, hits, hits_found,
                                     Comb_Approved, num_Comb,
                                     W_Comb);
            }
        }
    }

  return 0;
}

/*
  hits contains the indices that we want to estimate.

  x_samp contains a B-dimensional array for each of the `loops`
  iterations of the outer loop.  Every coordinate i of x "hashes to" a
  corresponding coordinate (permute[j] * i) mod B of x_samp[j], which
  gives an estimate of x[i].

  We estimate each coordinate as the median (independently in real and
  imaginary axes) of its images in the rows of x_samp.
 */
sfft_output
estimate_values(sfft_v1v2_data * data, const int *hits,
                const int &hits_found, complex_t * x_samp,
                const int &loops, int n, const int *permute, const int B,
                const int B2, const Filter & filter,
                const Filter & filter_Est, int location_loops)
{
  unsigned tid = omp_get_thread_num();
  sfft_v1v2_threadlocal_data *tldata = &data->threadlocal_data[tid];

  sfft_output ans;
  real_t **values = tldata->estimate_values_values;

  assert((unsigned)hits_found < data->x_samp_size);

  const __m128d signs = _mm_set_pd(-1.0, 1.0);

  for (int i = 0; i < hits_found; i++)
    {
      int position = 0;

      for (int j = 0; j < loops; j++)
        {
          int offset = MIN(j, location_loops) * B + MAX(0,
                       j -
                       location_loops)
                       * B2;

          int cur_B = (j < location_loops) ? B : B2;
          const Filter & cur_filter =
            (j < location_loops) ? filter : filter_Est;
          int permuted_index = timesmod(permute[j], hits[i], n);
          int hashed_to = permuted_index / (n / cur_B);
          int dist = permuted_index % (n / cur_B);
          if (dist > (n / cur_B) / 2)
            {
              hashed_to = (hashed_to + 1) % cur_B;
              dist -= n / cur_B;
            }
          dist = (n - dist) % n;

          double *d_filter_freq = (double *)cur_filter.freq;
          double *d_x_samp = (double *)x_samp;

          __m128d ab =
            _mm_load_pd(d_x_samp + 2 * (offset + hashed_to));
          __m128d cd = _mm_load_pd(d_filter_freq + 2 * dist);
          __m128d dc = _mm_shuffle_pd(cd, cd, 1);

          __m128d ac_bd = _mm_mul_pd(ab, cd);
          __m128d ad_bc = _mm_mul_pd(ab, dc);
          __m128d mad_bc = _mm_mul_pd(ad_bc, signs);

          __m128d acpbd_bcmad = _mm_hadd_pd(ac_bd, mad_bc);

          __m128d cd_squares = _mm_mul_pd(cd, cd);
          __m128d cd_squares_sum =
            _mm_hadd_pd(cd_squares, cd_squares);

          __m128d r = _mm_div_pd(acpbd_bcmad, cd_squares_sum);

          _mm_storel_pd(values[0] + position, r);
          _mm_storeh_pd(values[1] + position, r);

          position++;
        }

      int location = (loops - 1) / 2;

      for (int a = 0; a < 2; a++)
        std::nth_element(values[a], values[a] + location,
                         values[a] + position);
      real_t realv = values[0][location];
      real_t imagv = values[1][location];
      ans[hits[i]] = realv + I * imagv;
    }

  FLOPCOUNT_INCREMENT(hits_found * loops * 12 + hits_found * 2);

  return ans;
}

/*
  Outer loop of the algorithm.

  If we are performing the Comb heuristic, first we do so.

  Then, `loops` times:
    choose a random permutation
    run inner_loop_locate
    if in the first location_loops loops, also run inner_loop_filter

  at the end, `hits` contains the coordinates that appear at least
  loop_threshold of location_loops times.  We estimate the values at
  these coordinates as the median of the images x_samp[loops].

  Returns a map from coordinates to estimates.
 */
sfft_output
outer_loop(sfft_v1v2_data * data, complex_t * origx, int n,
           const Filter & filter, const Filter & filter_Est, int B2,
           int num, int B, int W_Comb, int Comb_loops, int loop_threshold,
           int location_loops, int loops)
{

  unsigned tid = omp_get_thread_num();
  sfft_v1v2_threadlocal_data *tldata = &data->threadlocal_data[tid];

  /* storage */
  int *permute = tldata->permute;
  int *permuteb = tldata->permuteb;
  int *score = tldata->score;
  int *Comb_Approved = tldata->Comb_Approved;
  int *J = tldata->J;
  complex_t *x_samp = tldata->x_samp;
  int *hits = tldata->hits;

  memset(x_samp, 0, data->x_samp_size * sizeof(*x_samp));
  memset(score, 0, n * sizeof(*score));

  PROFILING_START_SECTION("Random Numbers");
  /* random numbers */
  int *a = (int *)malloc(loops * sizeof(int));
  int *b = (int *)malloc(loops * sizeof(int));
  int *ai = (int *)malloc(loops * sizeof(int));

  for (int i = 0; i < loops; i++)
    {
      a[i] = 0;
      while (gcd(a[i], n) != 1)
        {
          a[i] = int (random() % n);
        }
      ai[i] = mod_inverse(a[i], n);
      b[i] = 0;
    }
  PROFILING_END_SECTION();

  int hits_found = 0;

  // BEGIN Comb
  int num_Comb = num;

  PROFILING_START_SECTION("Comb filter");
  if (WITH_COMB)
    {
      for (int i = 0; i < Comb_loops; i++)
        Comb_Filt(data, origx, n, num, W_Comb,
                  Comb_Approved + i * num);
    }
  PROFILING_END_SECTION();

  PROFILING_START_SECTION("Radix Sort");
  if (Comb_loops > 1)
    {
      radix_sort(Comb_Approved, Comb_loops * num);
      int Last = 0;
      for (int i = 1; i < Comb_loops * num; i++)
        {
          if (Comb_Approved[i] != Comb_Approved[Last])
            Comb_Approved[++Last] = Comb_Approved[i];
        }
      num_Comb = Last + 1;
    }
  PROFILING_END_SECTION();

  if (WITH_COMB)
    {
      hits_found = num_Comb * (n / W_Comb);
      for (int j = 0; j < n / W_Comb; j++)
        for (int i = 0; i < num_Comb; i++)
          hits[j * num_Comb + i] =
            j * W_Comb + Comb_Approved[i];
    }
  //END Comb

  //BEGIN INNER LOOPS
  for (int i = 0; i < loops; i++)
    {
      permute[i] = ai[i];
      permuteb[i] = b[i];
    }

  PROFILING_START_SECTION("inner_loop_locate");
  inner_loop_locate(data, origx, n, filter, filter_Est,
                    num, B, B2,
                    a, ai, b,
                    x_samp, J, loops, location_loops,
                    loop_threshold, score, hits, hits_found,
                    Comb_Approved, num_Comb, W_Comb);
  PROFILING_END_SECTION();

  //END INNER LOOPS

  //BEGIN ESTIMATION
  PROFILING_START_SECTION("Estimate Values");
  sfft_output ans =
    estimate_values(data, hits, hits_found, x_samp, loops, n, permute,
                    B, B2,
                    filter, filter_Est, location_loops);
  PROFILING_END_SECTION();

  return ans;
}
