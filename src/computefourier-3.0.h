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

#ifndef COMPUTEFOURIER_3_H
#define COMPUTEFOURIER_3_H

#include "sfft.h"
#include "fft.h"
#include "common.h"

#include <complex.h>
#include <map>
#include "fftw.h"

#define OPTIMIZE_FFTW 0
//#define  WITH_MANSOUR 0

//Comments located in the cc file.
void
alternate_fft(sfft_v3_data * data, sfft_output * ans, complex_t * origx,
              int n, int k, int W_Man, int Man_loops, int B1, int Wind1,
              int Gauss_loops1, complex_t * filtert1, complex_t * filterf1,
              int B2, int Wind2, int Gauss_loops2, complex_t * filtert2,
              complex_t * filterf2);

#endif
