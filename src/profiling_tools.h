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

#ifndef PROFILING_TOOLS_H
#define PROFILING_TOOLS_H

#if defined(USE_PROFILING) && defined(USE_FLOPCOUNTER)

#include <iostream>
#include "flopcount.h"

static double profiling_time_start;
static double profiling_time_end;
static sfft_flopcount_t profiling_flops_start;
static sfft_flopcount_t profiling_flops_end;
static const char *profiling_section_name;

#define PROFILING_START_SECTION(SECTION_NAME) {	\
    profiling_section_name = SECTION_NAME;		\
    profiling_time_start = get_time();			\
    profiling_flops_start = sfft_flopcount_value();	\
  }

#define PROFILING_END_SECTION() {					\
    profiling_time_end = get_time();					\
    profiling_flops_end = sfft_flopcount_value();			\
    double prof_time = profiling_time_end - profiling_time_start;	\
    sfft_flopcount_t prof_flops = profiling_flops_end-profiling_flops_start; \
    double prof_perf = prof_flops/prof_time*1e-9;			\
    std::cout << "[PROFILING] SECTION '" << profiling_section_name << "'\n"; \
    std::cout << "[PROFILING]   * TIME        : " << prof_time << " s\n"; \
    std::cout << "[PROFILING]   * FLOPS       : " << prof_flops << " Flop\n"; \
    std::cout << "[PROFILING]   * PERFORMANCE : " << prof_perf << " GFlop/s\n"; \
  }

#else

#define PROFILING_START_SECTION(...)
#define PROFILING_END_SECTION(...)

#endif

#endif
