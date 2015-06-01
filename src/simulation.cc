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

#include <cassert>
#include <iostream>
#include <unistd.h>
#include <cstdlib>

#include "sfft.h"
#include "common.h"
#include "fftw.h"
#include "simulation.h"

void usage(char *basename)
{
  std::cout << "Usage: " << basename << " [OPTIONS]\n";
  std::cout << "Options:\n";
  std::cout << "  -n num     Set the problem size to num\n";
  std::cout << "  -k num     Set the number of frequencies to num\n";
  std::cout << "  -r num     Perform num repetitions of the experiment\n";
  std::cout << "  -o         Use FFTW_MEASURE instead of FFTW_ESTIMATE\n";
  std::cout <<
            "  -v version Use specific sFFT version (valid arguments: 1, 2 or 3)\n";
  std::cout << "  -h         Print this help message\n";
}

void simulation::setup(int argc, char **argv)
{
  n = 16384;
  k = 50;
  repetitions = 1;
  int version = 1;
  int fftw_opt = FFTW_ESTIMATE;

  char ch;
  while ((ch = getopt(argc, argv, "hton:k:r:v:")) != EOF)
    {
      switch (ch)
        {
        case 'n':
          n = atoi(optarg);
          break;
        case 'k':
          k = atoi(optarg);
          break;
        case 'r':
          repetitions = atoi(optarg);
          break;
        case 'v':
          version = atoi(optarg);
          break;
        case 'o':
          fftw_opt = FFTW_MEASURE;
          break;
        case 'h':
        default:
          usage(argv[0]);
          exit(1);
        }
    }

  if (version == 1)
    this->plan = sfft_make_plan(n, k, SFFT_VERSION_1, fftw_opt);
  else if (version == 2)
    this->plan = sfft_make_plan(n, k, SFFT_VERSION_2, fftw_opt);
  else if (version == 3)
    this->plan = sfft_make_plan(n, k, SFFT_VERSION_3, fftw_opt);
  else
    {
      usage(argv[0]);
      exit(1);
    }

  generate_input();
}

void simulation::generate_input()
{
  // Create datastructures
  x = (complex_t *) sfft_malloc(n * sizeof(*x));

  srand(17);
  srand48(time(NULL) ^ (getpid() * 171717));

  //Randomized the None Zero Bins
  x_f = (complex_t *) calloc(n, sizeof(*x_f));
  int *LARGE_FREQ = (int *)malloc(k * sizeof(*LARGE_FREQ));
  for (int i = 0; i < k; i++)
    {
      LARGE_FREQ[i] = (unsigned)floor(drand48() * n);
      x_f[LARGE_FREQ[i]] = 1.0;	// Will ADD Random Phase and Amplitude Later.
    }
  fftw_dft(x, n, x_f, 1);
}

sfft_output simulation::run()
{
  sfft_output result;
  for (int i = 0; i < repetitions; i++)
    sfft_exec(this->plan, this->x, &result);

  return result;
}

complex_t *simulation::frequencies()
{
  return this->x_f;
}

int simulation::get_n()
{
  return this->n;
}

int simulation::get_k()
{
  return this->k;
}
