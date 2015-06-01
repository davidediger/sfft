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

#ifndef GAUSSFILTER_H
#define GAUSSFILTER_H

#include "fft.h"

struct Filter
{
  complex_t *time;
  int sizet;
  complex_t *freq;
};

/*
  Create a window function such that:
      the main lobe has width 2 * n * filterfrac
      outside the main lobe, the linf residual is tolerance
  Computes the required w.
  Allocates and returns the filter.
 */
complex_t *make_dolphchebyshev_t(double lobefrac, double tolerance, int &w);

complex_t *make_gaussian_t(double lobefrac, double tolerance, int &w);

complex_t *make_kaiserbessel_t(double lobefrac, double tolerance, int &w);

/*
  Modifies a w-dimensional window function to have n-dimensional FFT
  the sum of b adjacent ones previously.

  Allocates and returns a Filter instance pointing to the modified x and an n-dimensional FFT of it.
 */

Filter make_multiple_t(complex_t * x, int w, int n, int b);

#endif
