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

#ifndef UTILS_H
#define UTILS_H

#include "fft.h"

// Compute the gcd of a and b
// assumes a, b > 0.
int gcd(int a, int b);

// crappy inversion code I stole from elsewhere
// Undefined if gcd(a, n) > 1
int mod_inverse(int a, int n);

/*
  Compute the num'th smallest element of the length n input.
  uses std::nth_element, but doesn't mutate input.
 */
real_t nth_element_immutable(real_t * input, int n, int num, real_t * tmp);
int nth_int_element_immutable(int *input, int n, int num, int *tmp);

double mean(real_t * x, int n);

double variance(real_t * x, int n);

/*
  Output the indices corresponding to the num largest elements of samples.
  Output is sorted.
*/
void find_largest_indices(int *output, int num, real_t * samples, int n,
                          real_t * tmp_storage);

void radix(int byte, int size, int *A, int *TEMP);

void radix_sort(int *A, int size);

void radix_filt(int byte, int size, int *A, int *TEMP, complex_t * Filter,
                complex_t * TMP_F);

void radix_sort_filt(int *A, complex_t * Filter, int size);

int floor_to_pow2(double x);

//x[i] <- x[i-r]
void shift(complex_t * x, int n, int r);

double phase(complex_t x);

double AWGN(complex_t * x, int n, double std_noise);

inline double cabs2(complex_t x)
{
  return (creal(x) * creal(x) + cimag(x) * cimag(x));
}

double binomial_cdf(double prob, int n, int needed);

#endif
