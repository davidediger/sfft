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

#include <cstdlib>
#include <cassert>
#include <cmath>

#include <omp.h>

#include "intrinsics.h"
#include "sfft.h"
#include "common.h"
#include "utils.h"
#include "parameters.h"
#include "computefourier-1.0-2.0.h"
#include "computefourier-3.0.h"

/** Declarations **************************************************************/

static sfft_v1v2_data *sfft_v1v2_make_plan(unsigned int n, unsigned int k,
    int fftw_opt);
static sfft_v1_data *sfft_v1_make_plan(unsigned int n, unsigned int k,
                                       int fftw_opt);
static sfft_v2_data *sfft_v2_make_plan(unsigned int n, unsigned int k,
                                       int fftw_opt);
static sfft_v3_data *sfft_v3_make_plan(unsigned int n, unsigned int k,
                                       int fftw_opt);

static void sfft_v1_free_plan(sfft_v1_data * data);
static void sfft_v2_free_plan(sfft_v2_data * data);
static void sfft_v3_free_plan(sfft_v3_data * data);

static void sfft_v1(unsigned int n, unsigned int k, sfft_v1_data * data,
                    complex_t * in, sfft_output * out);
static void sfft_v2(unsigned int n, unsigned int k, sfft_v2_data * data,
                    complex_t * in, sfft_output * out);
static void sfft_v3(unsigned int n, unsigned int k, sfft_v3_data * data,
                    complex_t * in, sfft_output * out);

/** Public API ****************************************************************/

void *sfft_malloc(size_t s)
{
  return _mm_malloc(s, 16);
}

void sfft_free(void *p)
{
  if (p != NULL)
    _mm_free(p);
}

sfft_plan *sfft_make_plan(int n, int k, sfft_version version,
                          int fftw_optimization)
{
  void *data;

  switch (version)
    {
    case SFFT_VERSION_1:
      data = (void *)sfft_v1_make_plan(n, k, fftw_optimization);
      break;
    case SFFT_VERSION_2:
      data = (void *)sfft_v2_make_plan(n, k, fftw_optimization);
      break;
    case SFFT_VERSION_3:
      data = (void *)sfft_v3_make_plan(n, k, fftw_optimization);
      break;
    default:
      return NULL;
    }

  if (data == NULL)
    return NULL;

  sfft_plan *plan = (sfft_plan *) malloc(sizeof(sfft_plan));
  plan->version = version;
  plan->n = n;
  plan->k = k;
  plan->data = data;

  return plan;
}

void sfft_free_plan(sfft_plan * plan)
{
  switch (plan->version)
    {
    case SFFT_VERSION_1:
      sfft_v1_free_plan((sfft_v1_data *) plan->data);
      break;
    case SFFT_VERSION_2:
      sfft_v2_free_plan((sfft_v2_data *) plan->data);
      break;
    case SFFT_VERSION_3:
      sfft_v3_free_plan((sfft_v3_data *) plan->data);
      break;
    }
}

void sfft_exec(sfft_plan * plan, complex_t * in, sfft_output * out)
{
  switch (plan->version)
    {
    case SFFT_VERSION_1:
      sfft_v1(plan->n, plan->k, (sfft_v1_data *) plan->data, in, out);
      break;
    case SFFT_VERSION_2:
      sfft_v2(plan->n, plan->k, (sfft_v2_data *) plan->data, in, out);
      break;
    case SFFT_VERSION_3:
      sfft_v3(plan->n, plan->k, (sfft_v3_data *) plan->data, in, out);
      break;
    }
}

void
sfft_exec_many(sfft_plan * plan, int num, complex_t ** in, sfft_output * out)
{
  #pragma omp parallel for
  for (int i = 0; i < num; i++)
    {
      sfft_exec(plan, in[i], out + i);
    }
}

/** Static Functions **********************************************************/

static void sfft_v1_free_plan(sfft_v1_data * data)
{
  sfft_free(data->inner_loop_locate_loop_sizes);

  sfft_v1v2_threadlocal_data *tl = data->threadlocal_data;
  sfft_free(tl->permute);
  sfft_free(tl->permuteb);
  sfft_free(tl->x_samp);
  sfft_free(tl->score);
  sfft_free(tl->hits);
  sfft_free(tl->Comb_Approved);
  sfft_free(tl->J);
  sfft_free(tl->inner_loop_locate_x_sampt);
  sfft_free(tl->inner_loop_locate_samples);
  sfft_free(tl->estimate_values_values);
  sfft_free(tl->nth_element_storage);
  sfft_free(tl->nth_int_element_storage);

  fftw_destroy_plan(tl->fftw_plan_inner_loop_locate_location);
  fftw_destroy_plan(tl->fftw_plan_inner_loop_locate_estimation);
}

static void sfft_v2_free_plan(sfft_v2_data * data)
{
  sfft_v1_free_plan(data);

  sfft_v1v2_threadlocal_data *tl = data->threadlocal_data;
  sfft_free(tl->comb_filt_x_sampt);
  sfft_free(tl->comb_filt_samples);
  sfft_free(tl->inner_loop_filter_comb_permuted_approved);

  fftw_destroy_plan(tl->fftw_plan_comb_filt);
}

static void sfft_v3_free_plan(sfft_v3_data * data)
{
  sfft_free(data->filtert1);
  sfft_free(data->filterf1);
  sfft_free(data->filtert2);
  sfft_free(data->filterf2);
  sfft_free(data->fftw_plan_gauss_sizes);
  sfft_free(data->fftw_plan_gauss_perm_sizes);
  sfft_free(data->fftw_plan_man_sizes);

  sfft_v3_threadlocal_data *tl = data->threadlocal_data;
  sfft_free(tl->gauss_samples);
  sfft_free(tl->gauss_perm_samples);
  sfft_free(tl->man_samples);
  sfft_free(tl->est_freqs);
  sfft_free(tl->est_values);
  sfft_free(tl->nth_element_storage);
  sfft_free(tl->nth_int_element_storage);
  fftw_destroy_plan(tl->fftw_plan_gauss);
  fftw_destroy_plan(tl->fftw_plan_gauss_perm);
  fftw_destroy_plan(tl->fftw_plan_man);
}

static void
sfft_v1v2_create_threadlocal_data(unsigned n, unsigned k,
                                  sfft_v1v2_data * data,
                                  sfft_v1v2_threadlocal_data * tldata)
{
  int loops = data->loops_loc + data->loops_est;

  tldata->permute = (int *)sfft_malloc(loops * sizeof(int));
  tldata->permuteb = (int *)sfft_malloc(loops * sizeof(int));
  tldata->x_samp =
    (complex_t *) sfft_malloc(data->x_samp_size * sizeof(complex_t));
  tldata->score = (int *)sfft_malloc(n * sizeof(int));
  tldata->hits = (int *)sfft_malloc(n * sizeof(int));
  tldata->Comb_Approved =
    (int *)sfft_malloc(data->Comb_loops * data->B_thresh * sizeof(int));
  tldata->J = (int *)sfft_malloc(data->B_thresh * sizeof(int));
  tldata->inner_loop_locate_x_sampt =
    (complex_t *) sfft_malloc(data->x_samp_size * sizeof(complex_t));

  tldata->inner_loop_locate_samples =
    (real_t *) sfft_malloc(data->x_samp_size * sizeof(real_t));

  tldata->estimate_values_values =
    (real_t **) sfft_malloc(2 * sizeof(real_t *));
  for (unsigned a = 0; a < 2; a++)
    tldata->estimate_values_values[a] =
      (real_t *) sfft_malloc(loops * sizeof(real_t));

  if (WITH_COMB)
    {
      tldata->comb_filt_x_sampt =
        (complex_t *) sfft_malloc(data->W_Comb * sizeof(complex_t));
      tldata->comb_filt_samples =
        (real_t *) sfft_malloc(data->W_Comb * sizeof(real_t));
      int num_Comb = data->Comb_loops * data->B_thresh;	// we know num_Comb <= Comb_loops*num
      tldata->inner_loop_filter_comb_permuted_approved =
        (__typeof(tldata->inner_loop_filter_comb_permuted_approved))
        sfft_malloc(num_Comb *
                    sizeof
                    (*tldata->inner_loop_filter_comb_permuted_approved));
    }

  tldata->nth_element_storage =
    (real_t *) sfft_malloc(n * sizeof(real_t));
  tldata->nth_int_element_storage = (int *)sfft_malloc(n * sizeof(int));
}

static void
sfft_v1v2_plan_ffts(unsigned n, unsigned k, sfft_v1v2_data * data,
                    sfft_v1v2_threadlocal_data * tldata_vec)
{
  int threads = data->threads;
  int loc_loops = data->loops_loc;
  int est_loops = data->loops_est;

  if (WITH_COMB)
    {
      for (int i = 0; i < threads; i++)
        tldata_vec[i].fftw_plan_comb_filt =
          fftw_plan_dft_1d(data->W_Comb, (fftw_complex *)
                           tldata_vec[i].comb_filt_x_sampt,
                           (fftw_complex *)
                           tldata_vec[i].comb_filt_x_sampt,
                           FFTW_FORWARD, data->fftw_opt);
    }

  for (int i = 0; i < threads; i++)
    tldata_vec[i].fftw_plan_inner_loop_locate_location =
      fftw_plan_many_dft(1, data->inner_loop_locate_loop_sizes,
                         loc_loops, (fftw_complex *)
                         tldata_vec[i].inner_loop_locate_x_sampt,
                         NULL, 1, data->B_loc,
                         (fftw_complex *) tldata_vec[i].x_samp,
                         NULL, 1, data->B_loc, FFTW_FORWARD,
                         data->fftw_opt);

  for (int i = 0; i < threads; i++)
    tldata_vec[i].fftw_plan_inner_loop_locate_estimation =
      fftw_plan_many_dft(1,
                         data->inner_loop_locate_loop_sizes +
                         loc_loops, est_loops, (fftw_complex *)
                         tldata_vec[i].inner_loop_locate_x_sampt +
                         loc_loops * data->B_loc, NULL, 1,
                         data->B_est,
                         (fftw_complex *) tldata_vec[i].x_samp +
                         loc_loops * data->B_loc, NULL, 1,
                         data->B_est, FFTW_FORWARD,
                         data->fftw_opt);;
}

static sfft_v1v2_data *sfft_v1v2_make_plan(unsigned int n, unsigned int k,
    int fftw_opt)
{
  sfft_v1v2_data *data =
    (sfft_v1v2_data *) malloc(sizeof(sfft_v1v2_data));
  if (data == NULL)
    return NULL;

  double Bcst_loc = 1;
  double Bcst_est = 1;
  double Comb_cst = 2;
  int loc_loops = 4;
  int est_loops = 16;
  int threshold_loops = 3;
  data->Comb_loops = 1;
  double tolerance_loc = 1.e-8;
  double tolerance_est = 1.e-8;

  if (k > 50)
    get_expermient_vs_K_parameters(n, WITH_COMB, Bcst_loc, Bcst_est,
                                   Comb_cst, loc_loops, est_loops,
                                   threshold_loops,
                                   data->Comb_loops, tolerance_loc,
                                   tolerance_est);
  else
    get_expermient_vs_N_parameters(n, WITH_COMB, Bcst_loc, Bcst_est,
                                   Comb_cst, loc_loops, est_loops,
                                   threshold_loops,
                                   data->Comb_loops, tolerance_loc,
                                   tolerance_est);

  n = floor_to_pow2(n);
  assert(ALGORITHM1 || WITH_COMB);

  double BB_loc =
    (unsigned)(Bcst_loc * sqrt((double)(int)n * k / (log2(n))));
  double BB_est =
    (unsigned)(Bcst_est * sqrt((double)(int)n * k / (log2(n))));

  double lobefrac_loc = 0.5 / (BB_loc);
  double lobefrac_est = 0.5 / (BB_est);

  int b_loc = int (1.2 * 1.1 * ((real_t) n / BB_loc));
  int b_est = int (1.4 * 1.1 * ((real_t) n / BB_est));

  data->fftw_opt = fftw_opt;

  data->B_loc = floor_to_pow2(BB_loc);
  data->B_thresh = 2 * k;
  data->B_est = floor_to_pow2(BB_est);

  data->W_Comb = floor_to_pow2(Comb_cst * n / data->B_loc);

  data->loops_loc = loc_loops;
  data->loops_thresh = threshold_loops;
  data->loops_est = est_loops;

  // Create Filter
  int w_loc;
  complex_t *filtert =
    make_dolphchebyshev_t(lobefrac_loc, tolerance_loc, w_loc);
  data->filter = make_multiple_t(filtert, w_loc, n, b_loc);

  int w_est;
  complex_t *filtert_est =
    make_dolphchebyshev_t(lobefrac_est, tolerance_est, w_est);
  data->filter_est = make_multiple_t(filtert_est, w_est, n, b_est);

  // storage data
  int loops = data->loops_loc + data->loops_est;
  data->threads = omp_get_max_threads();
  data->fftw_opt = fftw_opt;

  data->x_samp_size =
    data->loops_loc * data->B_loc + data->loops_est * data->B_est;

  data->inner_loop_locate_loop_sizes =
    (int *)sfft_malloc(loops * sizeof(int));
  for (int i = 0; i < loc_loops; i++)
    data->inner_loop_locate_loop_sizes[i] = data->B_loc;
  for (int i = loc_loops; i < loops; i++)
    data->inner_loop_locate_loop_sizes[i] = data->B_est;

  data->threadlocal_data =
    (sfft_v1v2_threadlocal_data *) malloc(data->threads *
                                          sizeof
                                          (sfft_v1v2_threadlocal_data));
  for (unsigned i = 0; i < data->threads; i++)
    sfft_v1v2_create_threadlocal_data(n, k, data,
                                      &data->threadlocal_data[i]);

  sfft_v1v2_plan_ffts(n, k, data, data->threadlocal_data);

  return data;
}

static sfft_v1_data *sfft_v1_make_plan(unsigned int n, unsigned int k,
                                       int fftw_opt)
{
  ALGORITHM1 = true;
  WITH_COMB = false;

  sfft_v1_data *data =
    (sfft_v1_data *) sfft_v1v2_make_plan(n, k, fftw_opt);
  if (data == NULL)
    return NULL;

  return data;
}

static sfft_v2_data *sfft_v2_make_plan(unsigned int n, unsigned int k,
                                       int fftw_opt)
{
  ALGORITHM1 = true;
  WITH_COMB = true;

  sfft_v2_data *data =
    (sfft_v2_data *) sfft_v1v2_make_plan(n, k, fftw_opt);
  if (data == NULL)
    return NULL;

  return data;
}

static void
sfft_v3_plan_ffts(unsigned n, unsigned k, sfft_v3_data * data,
                  sfft_v3_threadlocal_data * threaddata_vec)
{
  /* Mansour Filter DFT */
  data->fftw_plan_man_sizes =
    (int *)sfft_malloc(data->Man_loops * sizeof(int));
  for (int i = 0; i < data->Man_loops; i++)
    data->fftw_plan_man_sizes[i] = data->W_Man;

  for (unsigned i = 0; i < data->threads; i++)
    {
      threaddata_vec[i].fftw_plan_man =
        fftw_plan_many_dft(1, data->fftw_plan_man_sizes,
                           data->Man_loops, (fftw_complex *)
                           threaddata_vec[i].man_samples, NULL, 2,
                           1, (fftw_complex *)
                           threaddata_vec[i].man_samples, NULL, 2,
                           1, FFTW_FORWARD, data->fftw_opt);
    }

  /* Gauss Filter DFT */
  data->fftw_plan_gauss_sizes =
    (int *)sfft_malloc(data->Gauss_loops * sizeof(int));
  for (int i = 0; i < data->Gauss_loops; i++)
    data->fftw_plan_gauss_sizes[i] = data->B_g1;

  for (unsigned i = 0; i < data->threads; i++)
    {
      threaddata_vec[i].fftw_plan_gauss =
        fftw_plan_many_dft(1, data->fftw_plan_gauss_sizes,
                           data->Gauss_loops, (fftw_complex *)
                           threaddata_vec[i].gauss_samples, NULL, 2,
                           1, (fftw_complex *)
                           threaddata_vec[i].gauss_samples, NULL, 2,
                           1, FFTW_FORWARD, data->fftw_opt);
    }

  /* 2nd Gauss Filter DFT */
  data->fftw_plan_gauss_perm_sizes =
    (int *)sfft_malloc(data->Gauss2_loops * sizeof(int));
  for (int i = 0; i < data->Gauss2_loops; i++)
    data->fftw_plan_gauss_perm_sizes[i] = data->B_g2;

  for (unsigned i = 0; i < data->threads; i++)
    {
      threaddata_vec[i].fftw_plan_gauss_perm =
        fftw_plan_many_dft(1, data->fftw_plan_gauss_perm_sizes,
                           data->Gauss2_loops, (fftw_complex *)
                           threaddata_vec[i].gauss_perm_samples,
                           NULL, 2, 1, (fftw_complex *)
                           threaddata_vec[i].gauss_perm_samples,
                           NULL, 2, 1, FFTW_FORWARD,
                           data->fftw_opt);
    }
}

static void
sfft_v3_create_threadlocal_data(unsigned int n, unsigned int k,
                                sfft_v3_data * data,
                                sfft_v3_threadlocal_data * threaddata)
{
  threaddata->man_samples =
    (complex_t *) sfft_malloc(data->Man_loops * data->W_Man *
                              sizeof(complex_t));
  threaddata->gauss_samples =
    (complex_t *) sfft_malloc(data->Gauss_loops * data->B_g1 *
                              sizeof(complex_t));
  threaddata->gauss_perm_samples =
    (complex_t *) sfft_malloc(data->Gauss2_loops * data->B_g2 *
                              sizeof(complex_t));

  threaddata->est_freqs = (int *)sfft_malloc(5 * k * sizeof(int));
  threaddata->est_values =
    (complex_t *) sfft_malloc(5 * k * sizeof(complex_t));
  threaddata->perm_x =
    (complex_t *) sfft_malloc(data->w_g2 * sizeof(complex_t));

  threaddata->nth_element_storage =
    (real_t *) sfft_malloc(n * sizeof(real_t));
  threaddata->nth_int_element_storage =
    (int *)sfft_malloc(n * sizeof(int));
}

static sfft_v3_data *sfft_v3_make_plan(unsigned int n, unsigned int k,
                                       int fftw_opt)
{
  sfft_v3_data *data = (sfft_v3_data *) sfft_malloc(sizeof(sfft_v3_data));
  if (data == NULL)
    return NULL;

  double Bcst = 1;
  double Bcst2 = 0.25;
  double Mancst = 10;

  data->Gauss_loops = 2;
  data->Gauss2_loops = 2;
  data->Man_loops = 2;

  int gwidth1 = 1;
  int gwidth2 = 1;
  double tolerance_g1 = 1.e-8;
  double tolerance_g2 = 1.e-8;

  n = floor_to_pow2(n);
  data->W_Man = floor_to_pow2(Mancst * ((double)k));
  data->W_Man = std::min(n / 2, (unsigned int)data->W_Man);

  //== First GAUSSIAN FILTER ===========================================

  real_t BB = (unsigned)(Bcst * ((double)k));
  data->B_g1 = floor_to_pow2(BB);

  double lobefrac_g1 = 0.5 / (BB);
  int b_g1 = gwidth1;

  if (b_g1 == 1)
    b_g1 = int (1.00 * ((double)n / data->B_g1));

  data->filtert1 =
    make_dolphchebyshev_t(lobefrac_g1, tolerance_g1, data->w_g1);
  data->filterf1 =
    make_multiple_t(data->filtert1, data->w_g1, n, b_g1).freq;

  //== Second GAUSSIAN FILTER ==========================================

  real_t BB2 = (unsigned)(Bcst2 * ((double)k));
  data->B_g2 = floor_to_pow2(BB2);

  double lobefrac_g2 = 0.5 / (BB2);
  int b_g2 = gwidth2;

  if (b_g2 == 1)
    b_g2 = int (1.00 * ((double)n / data->B_g2));

  data->filtert2 =
    make_dolphchebyshev_t(lobefrac_g2, tolerance_g2, data->w_g2);
  data->filterf2 =
    make_multiple_t(data->filtert2, data->w_g2, n, b_g2).freq;

  //== Storage =========================================================
  data->threads = omp_get_max_threads();
  data->fftw_opt = fftw_opt;

  data->threadlocal_data =
    (sfft_v3_threadlocal_data *) sfft_malloc(data->threads *
        sizeof
        (sfft_v3_threadlocal_data));
  for (unsigned i = 0; i < data->threads; i++)
    {
      sfft_v3_create_threadlocal_data(n, k, data,
                                      &data->threadlocal_data[i]);
    }

  sfft_v3_plan_ffts(n, k, data, data->threadlocal_data);

  return data;
}

static void
sfft_v1(unsigned int n, unsigned int k, sfft_v1_data * data,
        complex_t * in, sfft_output * out)
{
  *out = outer_loop(data, in, n, data->filter, data->filter_est,
                    data->B_est, data->B_thresh, data->B_loc,
                    data->W_Comb, data->Comb_loops, data->loops_thresh,
                    data->loops_loc, data->loops_loc + data->loops_est);
}

static void
sfft_v2(unsigned int n, unsigned int k, sfft_v2_data * data,
        complex_t * in, sfft_output * out)
{
  ALGORITHM1 = true;
  WITH_COMB = true;

  *out = outer_loop(data, in, n, data->filter, data->filter_est,
                    data->B_est, data->B_thresh, data->B_loc,
                    data->W_Comb, data->Comb_loops, data->loops_thresh,
                    data->loops_loc, data->loops_loc + data->loops_est);

}

static void
sfft_v3(unsigned int n, unsigned int k, sfft_v3_data * data,
        complex_t * in, sfft_output * out)
{
  alternate_fft(data, out, in, n, k, data->W_Man, data->Man_loops,
                data->B_g1, data->w_g1, data->Gauss_loops,
                data->filtert1, data->filterf1, data->B_g2, data->w_g2,
                data->Gauss_loops, data->filtert2, data->filterf2);

}
