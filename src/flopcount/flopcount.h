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

#ifndef SFFT_FLOPCOUNT_H
#define SFFT_FLOPCOUNT_H

#define FLOPCOUNT_COST_COMPLEX_ADD (2)
#define FLOPCOUNT_COST_COMPLEX_SUB (2)
#define FLOPCOUNT_COST_COMPLEX_MUL (6)
#define FLOPCOUNT_COST_COMPLEX_DIV (30)
#define FLOPCOUNT_COST_SQRT (1)
#define FLOPCOUNT_COST_SINCOS (29)
#define FLOPCOUNT_COST_PHASE (42)
#define FLOPCOUNT_COST_SINCOS_IPP (19)
#define FLOPCOUNT_COST_PHASE_IPP (37)

#ifdef USE_FLOPCOUNTER

/* The flopcount type */
typedef unsigned long long sfft_flopcount_t;

/* The actual instruction counter value */
extern sfft_flopcount_t _sfft_flopcount;

/* Returns the current value of the instruction counter */
sfft_flopcount_t sfft_flopcount_value();

#define FLOPCOUNT_RESET()   (_sfft_flopcount=0)
#define FLOPCOUNT_INCREMENT(n)  (_sfft_flopcount+=(n))
#define FLOPCOUNT_INCREMENT_FFTW_PLAN(p) {\
  double add, mul, fma;\
  fftw_flops(p, &add, &mul, &fma);\
  _sfft_flopcount += (int)add+mul+2*fma;\
}

#else

#define FLOPCOUNT_RESET(...)
#define FLOPCOUNT_INCREMENT(...)
#define FLOPCOUNT_INCREMENT_FFTW_PLAN(...)

#endif

#endif
