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

#ifndef FFT_H
#define FFT_H

#include<cmath>
#include<complex.h>
#include<vector>

#define Vec(a, b) std::vector<__typeof(*(a))> ((a), (a)+(b))

// allow easy change to float or long double
//#define USE_FLOAT
#define USE_DOUBLE

#ifdef USE_FLOAT
typedef float complex complex_t;
typedef float real_t;
#define cexp cexpf
#define exp expf
#endif

#ifdef USE_DOUBLE
typedef double complex complex_t;
typedef double real_t;
#endif

//#define DEBUG

#endif
