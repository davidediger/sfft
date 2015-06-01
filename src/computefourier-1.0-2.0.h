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

#ifndef COMPUTEFOURIER_1_2_H
#define COMPUTEFOURIER_1_2_H

#include "sfft.h"
#include "fft.h"
#include "common.h"

#include <complex.h>
#include "fftw.h"
#include "filters.h"

#define OPTIMIZE_FFTW 0

//Comments located in the cc file.
sfft_output
outer_loop(sfft_v1v2_data * data, complex_t * origx, int n,
           const Filter & filter, const Filter & filter_Est, int B2,
           int num, int B, int W_Comb, int Comb_loops, int loop_threshold,
           int location_loops, int loops);

#endif
