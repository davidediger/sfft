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

#ifndef RUN_SIMULATION_H
#define RUN_SIMULATION_H

#include "fft.h"
#include "sfft.h"

void usage(char *basename);

class simulation
{
protected:
  int n;
  int k;
  int repetitions;

  sfft_plan *plan;

  complex_t *x;
  complex_t *x_f;

public:
  simulation()
  {
  } simulation(int argc, char **argv)
  {
    setup(argc, argv);
  }

  complex_t *frequencies();
  int get_n();
  int get_k();

  void setup(int argc, char **argv);
  void generate_input();
  sfft_output run();
};

#endif
