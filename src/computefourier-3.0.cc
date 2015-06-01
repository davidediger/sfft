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

#ifdef HAVE_IPP
#include <ippvm.h>
#include <ipps.h>
#endif
#include <omp.h>

#include "intrinsics.h"
#include "computefourier-3.0.h"
#include "filters.h"
#include "utils.h"
#include "timer.h"

#include "flopcount.h"
#include "profiling_tools.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#ifdef HAVE_IPP
#define approx_atan2_vec2(a, b, c) ippsAtan2_64f_A26(a, b, c, 2)
#define approx_sincos(a, s, c) ippsSinCos_64f_A26 (a, s, c, 1)
#define approx_sincos_vec2(a, s, c) ippsSinCos_64f_A26 (a, s, c, 2)
#define approx_sqrt(s) ippsSqrt_64f_I(s, 1)
#define FLOPCOUNT_COST_APPROX_SINCOS (FLOPCOUNT_COST_SINCOS_IPP)
#define FLOPCOUNT_COST_APPROX_PHASE (FLOPCOUNT_COST_PHASE_IPP)
#else
#warning "Compiling without Intel Performance Primitives - resulting code will be slower!"
#define approx_atan2_vec2(a, b, c)  c[0] = atan2(a[0], b[0]); c[1] = atan2(a[1], b[1])
#define approx_sincos(a, s, c) *s = sin(*a); *c = cos(*a)
#define approx_sincos_vec2(a, s, c) s[0] = sin(a[0]); c[0] = cos(a[0]); s[1] = sin(a[1]); c[1] = cos(a[1])
#define approx_sqrt(s) *s = sqrt(*s)
#define FLOPCOUNT_COST_APPROX_SINCOS (FLOPCOUNT_COST_SINCOS)
#define FLOPCOUNT_COST_APPROX_PHASE (FLOPCOUNT_COST_PHASE)
#endif

// set dest such that in fourier space,
// F(dest)[i] has = F(src)[i ai + b]
// set F(x[i]) = F(x)[i ai + b] cyclically
static int
permute_fourier(complex_t * dest, int w, complex_t * src, int n,
                int offset, int ai, int b)
{
  complex_t shift =
    cos(2 * M_PI * offset * b / n) + I * sin(2 * M_PI * offset * b / n);
  complex_t step =
    cos(2 * M_PI * ai * b / n) + I * sin(2 * M_PI * ai * b / n);
  FLOPCOUNT_INCREMENT(2 *
                      (2 * FLOPCOUNT_COST_SINCOS + 9 +
                       FLOPCOUNT_COST_COMPLEX_MUL));
  int index = offset % n;

  for (int i = 0; i < w; i++)
    {
      dest[i] = src[index] * shift;
      shift = shift * step;
      FLOPCOUNT_INCREMENT(2 * FLOPCOUNT_COST_COMPLEX_MUL);
      index = (index + ai) % n;
    }
  return 0;
}

// Permute frequency
inline int timesmod(int x, int a, int n)
{
  return int ((((long long int)x) * a) % n);
}

#define Mansour_Filt Mansour_Filt_loops2
int
Mansour_Filt_loops2(sfft_v3_data * data, complex_t * origx, int n,
                    int W_Man, int init_offset, complex_t * x_man)
{

  assert(n % W_Man == 0);

  int n1 = n - 1;
  int sigma = n / W_Man;

  const unsigned tid = omp_get_thread_num();
  sfft_v3_threadlocal_data *tl_data = data->threadlocal_data + tid;

  complex_t *x_sampt = tl_data->man_samples;

  int index = init_offset;
  int index_p1 = (init_offset + 1) & n1;
  for (int i = 0; i < W_Man; i++)
    {
      __m128d o = _mm_load_pd((double *)(origx + index));
      __m128d p = _mm_load_pd((double *)(origx + index_p1));
      _mm_store_pd((double *)(x_sampt + 2 * i), o);
      _mm_store_pd((double *)(x_sampt + 2 * i + 1), p);
      index = (index + sigma) & n1;
      index_p1 = (index_p1 + sigma) & n1;
    }

  fftw_execute(tl_data->fftw_plan_man);
  FLOPCOUNT_INCREMENT_FFTW_PLAN(tl_data->fftw_plan_man);

  return 0;
}

#define Gauss_Filt Gauss_Filt_loops2
int
Gauss_Filt_loops2(sfft_v3_data * data, complex_t * origx, int n,
                  complex_t * filter, int w, int B, complex_t * x_gauss,
                  int init_G_offset)
{
  assert(n % B == 0);

  const int loops = data->Gauss_loops;

  const unsigned tid = omp_get_thread_num();
  sfft_v3_threadlocal_data *tl_data = data->threadlocal_data + tid;

  complex_t *x_sampt = tl_data->gauss_samples;
  memset(x_sampt, 0, loops * B * sizeof(complex_t));

  const double *d_origx = (double *)origx;
  const double *d_filter = (double *)filter;
  double *d_x_sampt = (double *)x_sampt;

  const unsigned n2_m_1 = 2 * n - 1;
  const unsigned origx_offset = (2 * init_G_offset + 2) & n2_m_1;

  const unsigned chunksize = 2 * B;
  const unsigned chunks = 2 * w / chunksize;

  for (unsigned chunk = 0; chunk < chunks; chunk++)
    {
      unsigned start = chunk * chunksize;
      unsigned end =
        std::min((chunk + 1) * chunksize, (unsigned)2 * w);

      __m128d a2b2 =
        _mm_load_pd(d_origx +
                    ((2 * init_G_offset + start) & n2_m_1));
      unsigned i2_mod_B = 0;
      for (unsigned i = start; i < end; i += 2)
        {
          __m128d ab = a2b2;
          a2b2 =
            _mm_load_pd(d_origx +
                        ((origx_offset + i) & n2_m_1));
          __m128d cd = _mm_load_pd(d_filter + i);

          __m128d cc = _mm_unpacklo_pd(cd, cd);
          __m128d dd = _mm_unpackhi_pd(cd, cd);

          __m128d a0a1 = _mm_unpacklo_pd(ab, a2b2);
          __m128d b0b1 = _mm_unpackhi_pd(ab, a2b2);

          __m128d ac = _mm_mul_pd(cc, a0a1);
          __m128d ad = _mm_mul_pd(dd, a0a1);
          __m128d bc = _mm_mul_pd(cc, b0b1);
          __m128d bd = _mm_mul_pd(dd, b0b1);

          __m128d ac_m_bd = _mm_sub_pd(ac, bd);
          __m128d ad_p_bc = _mm_add_pd(ad, bc);

          __m128d ab_times_cd = _mm_unpacklo_pd(ac_m_bd, ad_p_bc);
          __m128d a2b2_times_cd =
            _mm_unpackhi_pd(ac_m_bd, ad_p_bc);

          __m128d xy = _mm_load_pd(d_x_sampt + i2_mod_B);
          __m128d x2y2 = _mm_load_pd(d_x_sampt + i2_mod_B + 2);

          __m128d st = _mm_add_pd(xy, ab_times_cd);
          __m128d s2t2 = _mm_add_pd(x2y2, a2b2_times_cd);

          _mm_store_pd(d_x_sampt + i2_mod_B, st);
          _mm_store_pd(d_x_sampt + i2_mod_B + 2, s2t2);

          i2_mod_B += 4;
        }
    }

  FLOPCOUNT_INCREMENT(w * 2 *
                      (FLOPCOUNT_COST_COMPLEX_MUL +
                       FLOPCOUNT_COST_COMPLEX_ADD));

  fftw_execute(tl_data->fftw_plan_gauss);
  FLOPCOUNT_INCREMENT_FFTW_PLAN(tl_data->fftw_plan_gauss);

  return 0;
}

#define Gauss_Filt_Perm Gauss_Filt_Perm_loops2
int
Gauss_Filt_Perm_loops2(sfft_v3_data * data, complex_t * origx, int n,
                       complex_t * filter, int w, int B,
                       complex_t * x_gauss, int init_G_offset, int ai, int b0)
{
  assert(n % B == 0);
  assert(data->Gauss2_loops + w < n);

  const unsigned loops = 2;
  const unsigned tid = omp_get_thread_num();
  sfft_v3_threadlocal_data *tl_data = data->threadlocal_data + tid;

  complex_t *x_sampt = tl_data->gauss_perm_samples;
  complex_t *Perm_X = tl_data->perm_x;
  memset(x_sampt, 0, loops * B * sizeof(complex_t));

  assert((size_t) Perm_X % 16 == 0);
  assert((size_t) x_sampt % 16 == 0);

  int jump = 1;
  permute_fourier(Perm_X, w + jump * (loops - 1), origx, n,
                  init_G_offset, ai, b0);

  double *d_perm_x = (double *)Perm_X;
  double *d_filter = (double *)filter;
  double *d_x_sampt = (double *)x_sampt;

  const unsigned chunksize = 2 * B;
  const unsigned chunks = 2 * w / chunksize;

  for (unsigned chunk = 0; chunk < chunks; chunk++)
    {
      unsigned end =
        std::min((chunk + 1) * chunksize, (unsigned)2 * w);

      int i2_mod_B = 0;
      __m128d a2b2 = _mm_load_pd(d_perm_x + chunk * chunksize);
      for (unsigned i = chunk * chunksize; i < end; i += 2)
        {
          __m128d ab = a2b2;
          a2b2 = _mm_load_pd(d_perm_x + i + 2);
          __m128d cd = _mm_load_pd(d_filter + i);

          __m128d cc = _mm_unpacklo_pd(cd, cd);
          __m128d dd = _mm_unpackhi_pd(cd, cd);

          __m128d a0a1 = _mm_unpacklo_pd(ab, a2b2);
          __m128d b0b1 = _mm_unpackhi_pd(ab, a2b2);

          __m128d ac = _mm_mul_pd(cc, a0a1);
          __m128d ad = _mm_mul_pd(dd, a0a1);
          __m128d bc = _mm_mul_pd(cc, b0b1);
          __m128d bd = _mm_mul_pd(dd, b0b1);

          __m128d ac_m_bd = _mm_sub_pd(ac, bd);
          __m128d ad_p_bc = _mm_add_pd(ad, bc);

          __m128d ab_times_cd = _mm_unpacklo_pd(ac_m_bd, ad_p_bc);
          __m128d a2b2_times_cd =
            _mm_unpackhi_pd(ac_m_bd, ad_p_bc);

          __m128d xy = _mm_load_pd(d_x_sampt + i2_mod_B);
          __m128d x2y2 = _mm_load_pd(d_x_sampt + i2_mod_B + 2);

          __m128d st = _mm_add_pd(xy, ab_times_cd);
          __m128d s2t2 = _mm_add_pd(x2y2, a2b2_times_cd);

          _mm_store_pd(d_x_sampt + i2_mod_B, st);
          _mm_store_pd(d_x_sampt + i2_mod_B + 2, s2t2);

          i2_mod_B += 4;
        }
    }

  FLOPCOUNT_INCREMENT(w * loops *
                      (FLOPCOUNT_COST_COMPLEX_MUL +
                       FLOPCOUNT_COST_COMPLEX_ADD));

  fftw_execute(tl_data->fftw_plan_gauss_perm);
  FLOPCOUNT_INCREMENT_FFTW_PLAN(tl_data->fftw_plan_gauss_perm);

  return 0;
}

// REMOVE THE EFFECT OF A LARGE FREQUENCY FROM THE OUTPUT OF MANSOUR FILTER
static void
update_mansour_loops2(int key, complex_t value, complex_t * SAMP, int n,
                      int W_Man, int offset)
{
  const real_t PI2_OVER_N = 2 * M_PI / n;
  double *d_SAMP = (double *)SAMP;
  int hashed_to = key & (W_Man - 1);

  double value_real = creal(value);
  double value_imag = cimag(value);

  double cp_arg[] =
  { PI2_OVER_N * key * offset, PI2_OVER_N * key * (offset + 1) };
  double correct_phase_real[2];
  double correct_phase_imag[2];
  approx_sincos_vec2(cp_arg, correct_phase_imag, correct_phase_real);
  FLOPCOUNT_INCREMENT(4 * FLOPCOUNT_COST_APPROX_SINCOS + 4 /* MULTS */  +
                      1 /* ADDS */ );

  __m128d W_Man_vec = _mm_set1_pd(W_Man);
  __m128d value_real_vec = _mm_set1_pd(value_real);
  __m128d value_imag_vec = _mm_set1_pd(value_imag);
  __m128d cp_real_vec =
    _mm_set_pd(correct_phase_real[1], correct_phase_real[0]);
  __m128d cp_imag_vec =
    _mm_set_pd(correct_phase_imag[1], correct_phase_imag[0]);

  __m128d a0 = _mm_mul_pd(value_real_vec, cp_real_vec);
  __m128d b0 = _mm_mul_pd(value_imag_vec, cp_imag_vec);
  __m128d c0 = _mm_sub_pd(a0, b0);
  __m128d remove_real_vec = _mm_mul_pd(W_Man_vec, c0);
  __m128d a1 = _mm_mul_pd(value_real_vec, cp_imag_vec);
  __m128d b1 = _mm_mul_pd(value_imag_vec, cp_real_vec);
  __m128d c1 = _mm_add_pd(a1, b1);
  __m128d remove_imag_vec = _mm_mul_pd(W_Man_vec, c1);

  FLOPCOUNT_INCREMENT(4 * FLOPCOUNT_COST_COMPLEX_MUL);

  __m128d remove_0_vec =
    _mm_unpacklo_pd(remove_real_vec, remove_imag_vec);
  __m128d remove_1_vec =
    _mm_unpackhi_pd(remove_real_vec, remove_imag_vec);

  __m128d x0 = _mm_load_pd(d_SAMP + 4 * hashed_to);
  __m128d y0 = _mm_load_pd(d_SAMP + 4 * hashed_to + 2);
  __m128d x1 = _mm_sub_pd(x0, remove_0_vec);
  __m128d y1 = _mm_sub_pd(y0, remove_1_vec);
  _mm_store_pd(d_SAMP + 4 * hashed_to, x1);
  _mm_store_pd(d_SAMP + 4 * hashed_to + 2, y1);

  FLOPCOUNT_INCREMENT(2 * FLOPCOUNT_COST_COMPLEX_SUB);
}

// REMOVE THE EFFECT OF A LARGE FREQUENCY FROM THE OUTPUT OF GAUSSIAN FILTER
static void
update_gaussian_loops2(int key, complex_t value, complex_t * GAUSS_SAMP,
                       complex_t * filterf, int n, int B,
                       int init_G_offset, int step)
{
  const real_t PI2_DIV_N = 2 * M_PI / n;
  const int n1 = n - 1;
  FLOPCOUNT_INCREMENT(2);

  unsigned n_over_B = n / B;

  int hashed_to = int (key * 1. / n_over_B + 0.5) % B;
  FLOPCOUNT_INCREMENT(3);

  double *d_filterf = (double *)filterf;
  double *d_SAMP = (double *)GAUSS_SAMP;

  /* value */
  double a1 = ((double *)&value)[0];
  double b1 = ((double *)&value)[1];

  int dist1 = (hashed_to * n_over_B - key + n) & n1;
  /*  int dist2 = (((hashed_to+1)&B1)*n_over_B - key + n)&n1;
     int dist3 = (((hashed_to+B-1)&B1)*n_over_B - key + n)&n1; */
  int dist2 = (dist1 + n + n_over_B) & n1;
  int dist3 = (dist1 + n + n - n_over_B) & n1;

  int key_offset = (key * init_G_offset) % n;
  int key_offset2 = (key * (init_G_offset + step)) % n;

  /* correct phase */
  double cp_arg[] = { PI2_DIV_N * key_offset, PI2_DIV_N * key_offset2 };
  double a_vec[2];
  double b_vec[2];
  approx_sincos_vec2(cp_arg, b_vec, a_vec);
  double a2 = a_vec[0];
  double b2 = b_vec[0];
  double a22 = a_vec[1];
  double b22 = b_vec[1];
  FLOPCOUNT_INCREMENT(2 + 4 * FLOPCOUNT_COST_APPROX_SINCOS);

  /* value*correct_phase */
  double value_correct_phase_real = a1 * a2 - b1 * b2;
  double value_correct_phase_imag = a1 * b2 + b1 * a2;
  FLOPCOUNT_INCREMENT(4 + 2);

  /* value*correct_phase */
  double value_correct_phase_real2 = a1 * a22 - b1 * b22;
  double value_correct_phase_imag2 = a1 * b22 + b1 * a22;
  FLOPCOUNT_INCREMENT(4 + 2);

  __m128d v1r = _mm_set1_pd(value_correct_phase_real);
  __m128d v1i = _mm_set1_pd(value_correct_phase_imag);
  __m128d v2r = _mm_set1_pd(value_correct_phase_real2);
  __m128d v2i = _mm_set1_pd(value_correct_phase_imag2);

  __m128d ab31 = _mm_load_pd(d_filterf + 2 * dist1);
  __m128d ba31 = _mm_shuffle_pd(ab31, ab31, 1);
  __m128d ab32 = _mm_load_pd(d_filterf + 2 * dist2);
  __m128d ba32 = _mm_shuffle_pd(ab32, ab32, 1);
  __m128d ab33 = _mm_load_pd(d_filterf + 2 * dist3);
  __m128d ba33 = _mm_shuffle_pd(ab33, ab33, 1);

  __m128d t1r = _mm_mul_pd(v1r, ab31);
  __m128d t1i = _mm_mul_pd(v1i, ba31);
  __m128d remove1 = _mm_addsub_pd(t1r, t1i);

  __m128d t2r = _mm_mul_pd(v1r, ab32);
  __m128d t2i = _mm_mul_pd(v1i, ba32);
  __m128d remove2 = _mm_addsub_pd(t2r, t2i);

  __m128d t3r = _mm_mul_pd(v1r, ab33);
  __m128d t3i = _mm_mul_pd(v1i, ba33);
  __m128d remove3 = _mm_addsub_pd(t3r, t3i);

  __m128d t4r = _mm_mul_pd(v2r, ab31);
  __m128d t4i = _mm_mul_pd(v2i, ba31);
  __m128d remove4 = _mm_addsub_pd(t4r, t4i);

  __m128d t5r = _mm_mul_pd(v2r, ab32);
  __m128d t5i = _mm_mul_pd(v2i, ba32);
  __m128d remove5 = _mm_addsub_pd(t5r, t5i);

  __m128d t6r = _mm_mul_pd(v2r, ab33);
  __m128d t6i = _mm_mul_pd(v2i, ba33);
  __m128d remove6 = _mm_addsub_pd(t6r, t6i);

  FLOPCOUNT_INCREMENT(6 * (4 + 2));

  __m128d d11 = _mm_load_pd(d_SAMP + 4 * hashed_to);
  __m128d d12 = _mm_sub_pd(d11, remove1);
  _mm_store_pd(d_SAMP + 4 * hashed_to, d12);
  __m128d d21 = _mm_load_pd(d_SAMP + 4 * ((hashed_to + 1) % B));
  __m128d d22 = _mm_sub_pd(d21, remove2);
  _mm_store_pd(d_SAMP + 4 * ((hashed_to + 1) % B), d22);
  __m128d d31 = _mm_load_pd(d_SAMP + 4 * ((hashed_to + B - 1) % B));
  __m128d d32 = _mm_sub_pd(d31, remove3);
  _mm_store_pd(d_SAMP + 4 * ((hashed_to - 1 + B) % B), d32);
  __m128d d41 = _mm_load_pd(d_SAMP + 2 + 4 * hashed_to);
  __m128d d42 = _mm_sub_pd(d41, remove4);
  _mm_store_pd(d_SAMP + 2 + 4 * hashed_to, d42);
  __m128d d51 = _mm_load_pd(d_SAMP + 2 + 4 * ((hashed_to + 1) % B));
  __m128d d52 = _mm_sub_pd(d51, remove5);
  _mm_store_pd(d_SAMP + 2 + 4 * ((hashed_to + 1) % B), d52);
  __m128d d61 = _mm_load_pd(d_SAMP + 2 + 4 * ((hashed_to + B - 1) % B));
  __m128d d62 = _mm_sub_pd(d61, remove6);
  _mm_store_pd(d_SAMP + 2 + 4 * ((hashed_to + B - 1) % B), d62);

  FLOPCOUNT_INCREMENT(12);
}

inline void
UPDATE_ALL(int freq, complex_t median_value, int n, int init_offset,
           int W_Man, int Man_loops, complex_t * MAN_SAMP,
           int init_G_offset, int jump, int B1, complex_t * filterf1,
           complex_t * GAUSS_SAMP, int Gauss_loops1, int a, int ai,
           int shift, int B2, complex_t * filterf2,
           complex_t * GAUSS_SAMP_PERM, int Gauss_loops2)
{
  update_gaussian_loops2((timesmod(freq, ai, n) + shift) % n,
                         median_value, GAUSS_SAMP_PERM, filterf2, n, B2,
                         a * init_G_offset, jump);

  update_gaussian_loops2(freq, median_value, GAUSS_SAMP, filterf1, n, B1,
                         init_G_offset, 1);

  update_mansour_loops2(freq, median_value, MAN_SAMP, n, W_Man,
                        init_offset);
}

static int
estimate_freq_gauss_loops2(sfft_v3_data * data, int WHICH_FILTER,
                           int BUCKETS, complex_t * SAMP, int init_offset,
                           int LOOPS, int n, int a, int b, int jump,
                           complex_t * filterf, int *EST_FREQS,
                           complex_t * EST_VALUES)
{
  assert(LOOPS == 2);

  const unsigned n1 = n - 1;
  const unsigned B1 = BUCKETS - 1;
  const real_t PI2 = 2 * M_PI;
  const real_t N_OVER_PI2 = n / PI2;
  const real_t PI2_A_OFFSET_OVER_N = PI2 * a * init_offset / n;
  const real_t BUCKETS_OVER_N = (real_t) BUCKETS / n;
  const int N_OVER_BUCKETS = n / BUCKETS;
  const real_t collision_threshold = 1e-12;

  real_t separation = 1. / jump;

  real_t *d_SAMP = (real_t *) SAMP;

  complex_t median_value;
  int freq = 0;
  int found = 0;

  real_t ZERO_BUCK_CHECK_CUTOFF = 1e-8;

  double zero_buck_check[2];

  for (int i = 0; i < BUCKETS; i += 2)
    {
      __m128d a0b0 = _mm_load_pd(d_SAMP + 4 * i);
      __m128d a1b1 = _mm_load_pd(d_SAMP + 4 * i + 2);
      __m128d a2b2 = _mm_load_pd(d_SAMP + 4 * i + 4);
      __m128d a3b3 = _mm_load_pd(d_SAMP + 4 * i + 6);

      __m128d a0b0_sq = _mm_mul_pd(a0b0, a0b0);
      __m128d a1b1_sq = _mm_mul_pd(a1b1, a1b1);
      __m128d a2b2_sq = _mm_mul_pd(a2b2, a2b2);
      __m128d a3b3_sq = _mm_mul_pd(a3b3, a3b3);
      FLOPCOUNT_INCREMENT(8);

      __m128d c0c1 = _mm_hadd_pd(a0b0_sq, a1b1_sq);
      __m128d c2c3 = _mm_hadd_pd(a2b2_sq, a3b3_sq);
      FLOPCOUNT_INCREMENT(4);

      __m128d zbc = _mm_hadd_pd(c0c1, c2c3);
      FLOPCOUNT_INCREMENT(1);

      _mm_store_pd(zero_buck_check, zbc);

      for (unsigned j = 0; j < 2; j++)
        {
          if (zero_buck_check[j] > ZERO_BUCK_CHECK_CUTOFF)
            {
              real_t a0 = d_SAMP[4 * i + 4 * j];
              real_t b0 = d_SAMP[4 * i + 4 * j + 1];
              real_t a1 = d_SAMP[4 * i + 4 * j + 2];
              real_t b1 = d_SAMP[4 * i + 4 * j + 3];

              real_t c0 = (a0 * a0 + b0 * b0);
              real_t c1 = (a1 * a1 + b1 * b1);
              FLOPCOUNT_INCREMENT(6);

              real_t atan_real[] = { a0, a1 };
              real_t atan_imag[] = { b0, b1 };
              real_t atan_result[2];
              approx_atan2_vec2(atan_imag, atan_real,
                                atan_result);
              real_t d0 = (real_t) atan_result[0];
              real_t d1 = (real_t) atan_result[1];
              FLOPCOUNT_INCREMENT(2 *
                                  FLOPCOUNT_COST_APPROX_PHASE);

              real_t median_phase = 0;
              double slope = 0;
              unsigned hashed_to = 0;
              int dist = 0;

              complex_t filter_value;

              real_t median_phase_cos;
              real_t median_phase_sin;
              real_t error_b = c1 / c0 - 1;
              real_t error = error_b * error_b;
              FLOPCOUNT_INCREMENT(3);

              error /= n;
              FLOPCOUNT_INCREMENT(1);

              if ((error < collision_threshold)
                  && (c0 > 0.01))
                {
                  slope = (d1 - d0) * separation;
                  FLOPCOUNT_INCREMENT(2);

                  freq = lrint(N_OVER_PI2 * slope) + n;
                  freq = freq & n1;
                  FLOPCOUNT_INCREMENT(4);

                  hashed_to =
                    lrint(freq * BUCKETS_OVER_N) & B1;
                  FLOPCOUNT_INCREMENT(3);

                  if (hashed_to == i + j)
                    {
                      median_phase =
                        d0 -
                        PI2_A_OFFSET_OVER_N * freq;
                      FLOPCOUNT_INCREMENT(2);

                      approx_sqrt(&c0);
                      approx_sincos(&median_phase,
                                    &median_phase_sin,
                                    &median_phase_cos);
                      FLOPCOUNT_INCREMENT
                      (FLOPCOUNT_COST_SQRT +
                       2 *
                       FLOPCOUNT_COST_APPROX_SINCOS);

                      median_value =
                        median_phase_cos +
                        I * median_phase_sin;
                      median_value =
                        c0 * median_value;
                      FLOPCOUNT_INCREMENT(2);

                      dist =
                        (hashed_to *
                         N_OVER_BUCKETS - freq +
                         n) & n1;
                      filter_value = filterf[dist];
                      FLOPCOUNT_INCREMENT(4);

                      median_value =
                        median_value / filter_value;
                      FLOPCOUNT_INCREMENT
                      (FLOPCOUNT_COST_COMPLEX_DIV);

                      freq =
                        (timesmod(freq, a, n) - b +
                         n) & n1;
                      FLOPCOUNT_INCREMENT(4);

                      EST_FREQS[found] = freq;
                      EST_VALUES[found] =
                        median_value;
                      found++;
                    }
                }
            }
        }
    }

  return found;
}

static int
estimate_freq_mansour_loops2(sfft_v3_data * data, int BUCKETS,
                             complex_t * SAMP, int init_offset, int LOOPS,
                             int n, int a, int b, int jump,
                             complex_t * filterf, int *EST_FREQS,
                             complex_t * EST_VALUES)
{
  const real_t PI2 = 2 * M_PI;
  const real_t N_OVER_PI2 = n / PI2;
  const real_t PI2_OVER_N = PI2 / n;
  const unsigned W_Man = BUCKETS;
  const unsigned FREQ_MASK = ((n - 1) & ~(W_Man - 1));

  assert(LOOPS == 2);

  real_t *d_SAMP = (real_t *) SAMP;
  real_t *d_EST_VALUES = (real_t *) EST_VALUES;

  int freq1, freq2, freq3;
  int found = 0;

  real_t NORM = 1. / BUCKETS;
  real_t NORM2 = NORM * NORM;
  real_t ZERO_BUCK_CHECK_CUTOFF = 1e-8;

  double zero_buck_check[2];

  __m128d norm2vec = _mm_set1_pd(NORM2);

  for (int i = 0; i < BUCKETS; i += 2)
    {
      __m128d a0b0 = _mm_load_pd(d_SAMP + 4 * i);
      __m128d a1b1 = _mm_load_pd(d_SAMP + 4 * i + 2);
      __m128d a2b2 = _mm_load_pd(d_SAMP + 4 * i + 4);
      __m128d a3b3 = _mm_load_pd(d_SAMP + 4 * i + 6);

      __m128d a0b0_sq = _mm_mul_pd(a0b0, a0b0);
      __m128d a1b1_sq = _mm_mul_pd(a1b1, a1b1);
      __m128d a2b2_sq = _mm_mul_pd(a2b2, a2b2);
      __m128d a3b3_sq = _mm_mul_pd(a3b3, a3b3);
      FLOPCOUNT_INCREMENT(8);

      __m128d c0c1 = _mm_hadd_pd(a0b0_sq, a1b1_sq);
      __m128d c0c1_normed = _mm_mul_pd(c0c1, norm2vec);
      __m128d c2c3 = _mm_hadd_pd(a2b2_sq, a3b3_sq);
      __m128d c2c3_normed = _mm_mul_pd(c2c3, norm2vec);
      FLOPCOUNT_INCREMENT(8);

      __m128d zbc = _mm_hadd_pd(c0c1_normed, c2c3_normed);
      FLOPCOUNT_INCREMENT(1);

      _mm_store_pd(zero_buck_check, zbc);

      for (unsigned j = 0; j < 2; j++)
        {
          if (zero_buck_check[j] > ZERO_BUCK_CHECK_CUTOFF)
            {
              real_t a0 = d_SAMP[4 * i + 4 * j];
              real_t b0 = d_SAMP[4 * i + 4 * j + 1];
              real_t a1 = d_SAMP[4 * i + 4 * j + 2];
              real_t b1 = d_SAMP[4 * i + 4 * j + 3];

              real_t c0 = (a0 * a0 + b0 * b0) * NORM2;
              real_t c1 = (a1 * a1 + b1 * b1) * NORM2;
              FLOPCOUNT_INCREMENT(8);

              real_t atan_real[] = { a0 * NORM, a1 * NORM };
              real_t atan_imag[] = { b0 * NORM, b1 * NORM };
              real_t atan_result[2];

              approx_atan2_vec2(atan_imag, atan_real,
                                atan_result);
              real_t d0 = (real_t) atan_result[0];
              real_t d1 = (real_t) atan_result[1];
              FLOPCOUNT_INCREMENT(4 /* MULTS */  +
                                  2 *
                                  FLOPCOUNT_COST_APPROX_PHASE);

              real_t collision_threshold =
                1e-10, median_phase = 0;
              double slope = 0;

              real_t median_abs = c0;
              real_t median_abs_inv = 1. / median_abs;
              FLOPCOUNT_INCREMENT(1);

              real_t b = (c1 * median_abs_inv) - 1;
              real_t error = b * b;
              FLOPCOUNT_INCREMENT(3);

              if ((error < n * collision_threshold)
                  && (median_abs > 0.01))
                {
                  slope = d1 - d0;
                  FLOPCOUNT_INCREMENT(1);

                  freq1 = lrint(slope * N_OVER_PI2);
                  freq2 = freq1 & FREQ_MASK;
                  freq3 = freq2 | (i + j);
                  FLOPCOUNT_INCREMENT(3);

                  int freq_offset = freq3 * init_offset;
                  median_phase =
                    d0 - PI2_OVER_N * freq_offset;
                  FLOPCOUNT_INCREMENT(4);

                  approx_sqrt(&median_abs);
                  real_t median_phase_cos;
                  real_t median_phase_sin;
                  approx_sincos(&median_phase,
                                &median_phase_sin,
                                &median_phase_cos);
                  real_t median_value_real =
                    median_abs * median_phase_cos;
                  real_t median_value_imag =
                    median_abs * median_phase_sin;
                  FLOPCOUNT_INCREMENT(FLOPCOUNT_COST_SQRT
                                      +
                                      2 *
                                      FLOPCOUNT_COST_APPROX_SINCOS
                                      + 2);

                  EST_FREQS[found] = freq3;
                  d_EST_VALUES[2 * found] =
                    median_value_real;
                  d_EST_VALUES[2 * found + 1] =
                    median_value_imag;
                  found++;
                }
            }
        }
    }

  return found;
}

/*
  Alternating FFT

  Return a map from coordinates to estimates.
*/

void
alternate_fft(sfft_v3_data * data, sfft_output * out, complex_t * origx,
              int n, int k, int W_Man, int Man_loops, int B1, int Wind1,
              int Gauss_loops1, complex_t * filtert1, complex_t * filterf1,
              int B2, int Wind2, int Gauss_loops2, complex_t * filtert2,
              complex_t * filterf2)
{
  const unsigned tid = omp_get_thread_num();
  sfft_v3_threadlocal_data *tl_data = data->threadlocal_data + tid;

  complex_t *MAN_SAMP = tl_data->man_samples;
  complex_t *GAUSS_SAMP = tl_data->gauss_samples;
  complex_t *GAUSS_SAMP_PERM = tl_data->gauss_perm_samples;

  sfft_output & ans = *out;

  int a = 0;
  while (gcd(a, n) != 1)
    {
      a = int (random() % n);
    }
  int ai = mod_inverse(a, n);
  int b = int (random() % n);
  int shift = ((ai * b) % n + n) % n;

  int init_offset = (unsigned)floor(drand48() * n);
  int init_G_offset = (unsigned)floor(drand48() * n);

//=======================================================================================================================================================//

  // #################################################################
  // ################                        #########################
  // ################    DO COMB FILTER      #########################
  // ################                        #########################
  // #################################################################

  int offset;

  PROFILING_START_SECTION("Do Comb Filter");
  Mansour_Filt(data, origx, n, W_Man, init_offset, MAN_SAMP);
  PROFILING_END_SECTION();

//=====================================================================================================================================================

  // #################################################################
  // ################                        #########################
  // ################  ESTIMATE FREQUENCIES  #########################
  // ################      COMB FILTER       #########################
  // #################################################################

  int found = 0;
  int found_now = 0;

  int *EST_FREQS = tl_data->est_freqs;
  complex_t *EST_VALUES = tl_data->est_values;

  PROFILING_START_SECTION("Estimate Frequencies of Comb Filter");

  found_now =
    estimate_freq_mansour_loops2(data, W_Man, MAN_SAMP, init_offset,
                                 Man_loops, n, 1, 0, 1, filterf1,
                                 EST_FREQS, EST_VALUES);

  for (int i = 0; i < found_now; i++)
    {

      ans[EST_FREQS[i]] = EST_VALUES[i];
      found++;

      for (int j = 0; j < Man_loops; j++)
        MAN_SAMP[j + 2 * (EST_FREQS[i] % W_Man)] = 0;
    }

  PROFILING_END_SECTION();

  if (found == k)
    return;

//=====================================================================================================================================================

  // #################################################################
  // ################                        #########################
  // ################  DO GAUSSIAN FILTER 1  #########################
  // ################                        #########################
  // #################################################################

  PROFILING_START_SECTION("Do Gaussian Filter 1");
  Gauss_Filt(data, origx, n, filtert1, Wind1, B1, GAUSS_SAMP,
             init_G_offset);
  PROFILING_END_SECTION();

  offset = init_G_offset;
  int key;
  complex_t value;
  int jump = 1;

  PROFILING_START_SECTION("UPDATE_GAUSSIAN");
  for (__typeof(ans.begin())it = ans.begin(); it != ans.end(); it++)
    {
      key = it->first;
      value = it->second;
      update_gaussian_loops2(key, value, GAUSS_SAMP, filterf1, n, B1,
                             offset, 1);
    }
  PROFILING_END_SECTION();

//=====================================================================================================================================================

  // #################################################################
  // ################                        #########################
  // ################ ESTIMATE FREQUENCIES   #########################
  // ################      GAUSSIAN 1        #########################
  // #################################################################

  PROFILING_START_SECTION("Estimate Frequencies of Gaussian Filter 1");
  found_now =
    estimate_freq_gauss_loops2(data, 0, B1, GAUSS_SAMP, init_G_offset,
                               Gauss_loops1, n, 1, 0, jump, filterf1,
                               EST_FREQS, EST_VALUES);
  PROFILING_END_SECTION();

  PROFILING_START_SECTION("Update Mansour & Gauss1");
  for (int i = 0; i < found_now; i++)
    {

      if (cabs2(ans[EST_FREQS[i]]) == 0)
        found++;
      FLOPCOUNT_INCREMENT(2 + 1);

      ans[EST_FREQS[i]] = ans[EST_FREQS[i]] + EST_VALUES[i];

      update_gaussian_loops2(EST_FREQS[i], EST_VALUES[i], GAUSS_SAMP,
                             filterf1, n, B1, init_G_offset, 1);
      update_mansour_loops2(EST_FREQS[i], EST_VALUES[i], MAN_SAMP, n,
                            W_Man, init_offset);
      /*
         UPDATE_ALL(EST_FREQS[i], EST_VALUES[i], n, init_offset, W_Man, Man_loops, MAN_SAMP,
         init_G_offset, jump,  B1, filterf1, GAUSS_SAMP, Gauss_loops1,
         a, ai, shift, B2, filterf2, GAUSS_SAMP_PERM, Gauss_loops2);
       */
    }
  PROFILING_END_SECTION();
//=====================================================================================================================================================

  // #################################################################
  // ################                        #########################
  // ################  DO GAUSSIAN FILTER 2  #########################
  // ################                        #########################
  // #################################################################

  PROFILING_START_SECTION("Do Gaussian Filter 2");
  Gauss_Filt_Perm(data, origx, n, filtert2, Wind2, B2, GAUSS_SAMP_PERM,
                  init_G_offset, ai, b);
  PROFILING_END_SECTION();

  PROFILING_START_SECTION("UPDATE_GAUSSIAN");
  for (__typeof(ans.begin())it = ans.begin(); it != ans.end(); it++)
    {
      key = it->first;
      value = it->second;
      update_gaussian_loops2((timesmod(key, ai, n) + shift) % n,
                             value, GAUSS_SAMP_PERM, filterf2, n, B2,
                             a * init_G_offset, jump);
    }
  PROFILING_END_SECTION();

//=====================================================================================================================================================

  // #################################################################
  // ################                        #########################
  // ################  ESTIMATE FREQUENCIES  #########################
  // ################      GAUSSIAN 2        #########################
  // ################                        #########################
  // #################################################################

  PROFILING_START_SECTION("Estimate Frequencies of Gaussian Filter 2");
  found_now =
    estimate_freq_gauss_loops2(data, 0, B2, GAUSS_SAMP_PERM,
                               init_G_offset, Gauss_loops2, n, a, b,
                               jump, filterf2, EST_FREQS, EST_VALUES);
  PROFILING_END_SECTION();

  PROFILING_START_SECTION("UPDATE_ALL");
  for (int i = 0; i < found_now; i++)
    {
      if (cabs2(ans[EST_FREQS[i]]) == 0)
        found++;
      FLOPCOUNT_INCREMENT(2 + 1);

      ans[EST_FREQS[i]] = ans[EST_FREQS[i]] + EST_VALUES[i];
      FLOPCOUNT_INCREMENT(FLOPCOUNT_COST_COMPLEX_ADD);

      UPDATE_ALL(EST_FREQS[i], EST_VALUES[i], n, init_offset, W_Man,
                 Man_loops, MAN_SAMP, init_G_offset, jump, B1,
                 filterf1, GAUSS_SAMP, Gauss_loops1, a, ai, shift, B2,
                 filterf2, GAUSS_SAMP_PERM, Gauss_loops2);

    }
  PROFILING_END_SECTION();

  //=====================================================================================================================================================
  // #################################################################
  // ################                        #########################
  // ################ LOOP BETWEEN FILTERS   #########################
  // ################                        #########################
  // #################################################################

  int nana = 0;

  int prev_count_man = 0, prev_count_gauss1 = 0, prev_count_gauss2 = 0;
  int count_man = 0, count_gauss1 = 0, count_gauss2 = 0;

  PROFILING_START_SECTION("Loop Between Filters");
  while (true)
    {

      if (nana % 3 == 0)
        found_now =
          estimate_freq_mansour_loops2(data, W_Man, MAN_SAMP,
                                       init_offset, Man_loops,
                                       n, 1, 0, 1, filterf1,
                                       EST_FREQS, EST_VALUES);
      else if (nana % 3 == 1)
        found_now =
          estimate_freq_gauss_loops2(data, 0, B1, GAUSS_SAMP,
                                     init_G_offset,
                                     Gauss_loops1, n, 1, 0,
                                     jump, filterf1,
                                     EST_FREQS, EST_VALUES);
      else
        found_now =
          estimate_freq_gauss_loops2(data, 0, B2,
                                     GAUSS_SAMP_PERM,
                                     init_G_offset,
                                     Gauss_loops2, n, a, b,
                                     jump, filterf2,
                                     EST_FREQS, EST_VALUES);

      for (int i = 0; i < found_now; i++)
        {

          if (cabs2(ans[EST_FREQS[i]]) == 0)
            found++;
          FLOPCOUNT_INCREMENT(2 + 1);

          ans[EST_FREQS[i]] = ans[EST_FREQS[i]] + EST_VALUES[i];
          FLOPCOUNT_INCREMENT(FLOPCOUNT_COST_COMPLEX_ADD);

          UPDATE_ALL(EST_FREQS[i], EST_VALUES[i], n, init_offset,
                     W_Man, Man_loops, MAN_SAMP, init_G_offset,
                     jump, B1, filterf1, GAUSS_SAMP, Gauss_loops1,
                     a, ai, shift, B2, filterf2, GAUSS_SAMP_PERM,
                     Gauss_loops2);

        }

      if (nana % 3 == 2)
        {

          count_man = 0;
          count_gauss1 = 0;
          count_gauss2 = 0;

          for (int j = 0; j < B1; j++)
            {
              if (cabs2(GAUSS_SAMP[j]) > 1e-6)
                count_gauss1++;
              FLOPCOUNT_INCREMENT(2 + 1);
            }
          for (int j = 0; j < B2; j++)
            {
              if (cabs2(GAUSS_SAMP_PERM[j]) > 1e-6)
                count_gauss2++;
              FLOPCOUNT_INCREMENT(2 + 1);
            }
          for (int j = 0; j < W_Man; j++)
            {
              if (cabs2(MAN_SAMP[2 * j] / W_Man) > 1e-6)
                count_man++;
              FLOPCOUNT_INCREMENT(2 + 1);
            }

          if (prev_count_man == count_man
              && prev_count_gauss1 == count_gauss1
              && prev_count_gauss2 == count_gauss2)
            break;
          prev_count_man = count_man;
          prev_count_gauss1 = count_gauss1;
          prev_count_gauss2 = count_gauss2;
        }

      nana++;
    }

  PROFILING_END_SECTION();

  return;
}
