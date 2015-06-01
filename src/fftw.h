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

#ifndef FFTW_H
#define FFTW_H

#include "fft.h"
#include <fftw3.h>

#ifdef USE_FLOAT

#define fftw_plan_dft_1d fftwf_plan_dft_1d
#define fftw_plan fftwf_plan
#define fftw_execute fftwf_execute
#define fftw_destroy_plan fftwf_destroy_plan
#define fftw_free fftwf_free

#endif

#define OPTIMIZE_FFTW 0
int fftw_dft(complex_t * out, int n, complex_t * x, int backwards = 0);

#endif
