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

#include "fft.h"
#include <algorithm>
#include <string.h>
#include <cassert>

//x[i] <- x[i-r]
void shift(complex_t * x, int n, int r)
{
  r = (n + r) % n;
  assert(n >= r);
  complex_t *tmp = (complex_t *) malloc(r * sizeof(*tmp));
  memcpy(tmp, x + n - r, r * sizeof(*tmp));
  memmove(x + r, x, (n - r) * sizeof(*x));
  memcpy(x, tmp, r * sizeof(*tmp));
  free(tmp);
}

// Compute the gcd of a and b
// assumes a, b > 0.
int gcd(int a, int b)
{
  if (a % b == 0)
    return b;
  return gcd(b, a % b);
}

double phase(complex_t x)
{
  return atan2(cimag(x), creal(x));
}

double mean(real_t * x, int n)
{

  real_t SUM = 0;

  for (int i = 0; i < n; i++)
    SUM += x[i];

  SUM /= n;

  return SUM;

}

double variance(real_t * x, int n)
{

  real_t avg = mean(x, n);
  real_t var = 0;

  for (int i = 0; i < n; i++)
    var += x[i] * x[i];

  var /= n;

  var -= avg * avg;

  return var;
}

// crappy inversion code I stole from elsewhere
// Undefined if gcd(a, n) > 1
int mod_inverse(int a, int n)
{
  int i = n, v = 0, d = 1;
  while (a > 0)
    {
      int t = i / a, x = a;
      a = i % x;
      i = x;
      x = d;
      d = v - t * x;
      v = x;
    }
  v %= n;
  if (v < 0)
    v = (v + n) % n;
  return v;
}

/*
  Compute the num'th smallest element of the length n input.
  uses std::nth_element, but doesn't mutate input.
*/
real_t
nth_element_immutable(real_t * input, int n, int num, real_t * tmp_storage)
{
  real_t *x = tmp_storage;
  memcpy(x, input, n * sizeof(*x));
  std::nth_element(x, x + num, x + n);
  real_t ans = x[num];
  return ans;
}

int nth_int_element_immutable(int *input, int n, int num, int *tmp_storage)
{
  int *x = tmp_storage;
  memcpy(x, input, n * sizeof(*x));
  std::nth_element(x, x + num, x + n);
  int ans = x[num];
  return ans;
}

/*
  Output the indices corresponding to the num largest elements of samples.
  Output is sorted.
*/
void
find_largest_indices(int *output, int num, real_t * samples, int n,
                     real_t * tmp_storage)
{
  assert(n >= num + 1);
  //use num+1 so we can use > cutoff and probably get exactly num.
  //if we get fewer, the second pass uses == cutoff.
  real_t cutoff =
    nth_element_immutable(samples, n, n - num - 1, tmp_storage);

  int count = 0;
  for (int i = 0; i < n; i++)
    if (samples[i] > cutoff)
      output[count++] = i;
  if (count < num)
    {
      for (int i = 0; i < n; i++)
        {
          if (samples[i] == cutoff)
            {
              output[count++] = i;
              if (count >= num)
                break;
            }
        }
      std::sort(output, output + count);
    }
  assert(count == num);
}

void radix(int byte, int size, int *A, int *TEMP)
{

  int *COUNT = (int *)calloc(256, sizeof(*COUNT));

  byte = byte << 3;

  for (int i = 0; i < size; ++i)
    ++COUNT[((A[i]) >> (byte)) & 0xFF];
  for (int i = 1; i < 256; ++i)
    COUNT[i] += COUNT[i - 1];
  for (int i = size - 1; i >= 0; --i)
    {
      TEMP[COUNT[(A[i] >> (byte)) & 0xFF] - 1] = A[i];
      --COUNT[(A[i] >> (byte)) & 0xFF];
    }

  free(COUNT);

}

void radix_sort(int *A, int size)
{
  int *TEMP = (int *)malloc(size * sizeof(*TEMP));

  for (unsigned int i = 0; i < sizeof(int); i += 2)
    {

      // even byte
      radix(i, size, A, TEMP);

      // odd byte
      radix(i + 1, size, TEMP, A);
    }

  free(TEMP);
}

void
radix_filt(int byte, int size, int *A, int *TEMP, complex_t * Filter,
           complex_t * TMP_F)
{

  int *COUNT = (int *)calloc(256, sizeof(*COUNT));

  byte = byte << 3;

  for (int i = 0; i < size; ++i)
    ++COUNT[((A[i]) >> (byte)) & 0xFF];
  for (int i = 1; i < 256; ++i)
    COUNT[i] += COUNT[i - 1];
  for (int i = size - 1; i >= 0; --i)
    {
      TEMP[COUNT[(A[i] >> (byte)) & 0xFF] - 1] = A[i];
      TMP_F[COUNT[(A[i] >> (byte)) & 0xFF] - 1] = Filter[i];
      --COUNT[(A[i] >> (byte)) & 0xFF];
    }

  free(COUNT);

}

void radix_sort_filt(int *A, complex_t * Filter, int size)
{

  int *TEMP = (int *)malloc(size * sizeof(*TEMP));

  complex_t *TMP_F = (complex_t *) malloc(size * sizeof(*TMP_F));

  for (unsigned int i = 0; i < sizeof(int); i += 2)
    {

      // even byte
      radix_filt(i, size, A, TEMP, Filter, TMP_F);

      // odd byte
      radix_filt(i + 1, size, TEMP, A, TMP_F, Filter);
    }

  free(TEMP);
  free(TMP_F);
}

int floor_to_pow2(double x)
{
  unsigned int ans;
  for (ans = 1; ans <= x; ans <<= 1) ;
  return ans / 2;
}

double AWGN(complex_t * x, int n, double std_noise)
{

  if (std_noise == 0)
    return 1000000000;

  complex_t gn = 0;

  double sig_power = 0;
  double noise_power = 0;
  double snr;
  double u, v;

  for (int h = 0; h < n; h++)
    {
      sig_power += cabs(x[h]) * cabs(x[h]);

      u = drand48();
      v = drand48();
      gn = std_noise * sqrt(-2 * log(u)) * cexp(2 * M_PI * I * v);

      noise_power += -2 * log(u);

      x[h] += gn;
    }

  noise_power = noise_power * std_noise * std_noise;
  snr = sig_power / noise_power;

  return snr;
}

double binomial_cdf(double prob, int n, int needed)
{
  double ans = 0;
  double choose = 1;
  for (int i = n; i >= needed; i--)
    {
      ans += choose * pow(prob, i) * pow(1 - prob, n - i);
      choose = choose * i / (n - i + 1);
    }
  return ans;
}
