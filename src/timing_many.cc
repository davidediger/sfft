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

#include <iostream>
#include <cassert>
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

#include <valgrind/callgrind.h>

#include "sfft.h"
#include "common.h"
#include "fftw.h"

static void usage(char *basename)
{
  std::cout << "Usage: " << basename << " [OPTIONS]\n";
  std::cout << "Options:\n";
  std::cout << "  -n num     Set the problem size to num\n";
  std::cout << "  -k num     Set the number of frequencies to num\n";
  std::cout << "  -s         Simple parallelism (no data sharing)\n";
  std::cout << "  -i num     Perform the sFFT on num inputs\n";
  std::cout << "  -o         Use FFTW_MEASURE instead of FFTW_ESTIMATE\n";
  std::cout <<
            "  -v version Use specific sFFT version (valid arguments: 1, 2 or 3)\n";
  std::cout << "  -h         Print this help message\n";
}

static void
parse_arguments(int argc, char **argv, int *n, int *k, int *num_inputs,
                int *version, int *fftw_opt, bool * simple_parallelism)
{
  char ch;
  while ((ch = getopt(argc, argv, "htosi:n:k:v:")) != EOF)
    {
      switch (ch)
        {
        case 'n':
          *n = atoi(optarg);
          break;
        case 'k':
          *k = atoi(optarg);
          break;
        case 's':
          *simple_parallelism = true;
          break;
        case 'i':
          *num_inputs = atoi(optarg);
          break;
        case 'v':
          *version = atoi(optarg);
          break;
        case 'o':
          *fftw_opt = FFTW_MEASURE;
          break;
        case 'h':
        default:
          usage(argv[0]);
          exit(1);
        }
    }
}

static sfft_plan *plan_by_params(int n, int k, int version, int fftw_opt)
{
  sfft_plan *plan = NULL;
  if (version == 1)
    plan = sfft_make_plan(n, k, SFFT_VERSION_1, fftw_opt);
  else if (version == 2)
    plan = sfft_make_plan(n, k, SFFT_VERSION_2, fftw_opt);
  else if (version == 3)
    plan = sfft_make_plan(n, k, SFFT_VERSION_3, fftw_opt);

  return plan;
}

static complex_t **generate_input(int n, int k, int num_inputs)
{
  srand(17);
  srand48(time(NULL) ^ (getpid() * 171717));

  complex_t **x = (complex_t **) malloc(num_inputs * sizeof(complex_t *));
  complex_t *x_f = (complex_t *) sfft_malloc(n * sizeof(*x_f));

  for (int t = 0; t < num_inputs; t++)
    {
      x[t] = (complex_t *) sfft_malloc(n * sizeof(complex_t));
      memset(x_f, 0, n * sizeof(complex_t));

      for (int i = 0; i < k; i++)
        {
          unsigned f = (unsigned)floor(drand48() * n);
          x_f[f] = 1.0;
        }
      fftw_dft(x[t], n, x_f, 1);
    }

  free(x_f);
  return x;
}

int main(int argc, char **argv)
{
  int n = 1 << 18;
  int k = 1 << 14;
  int num_inputs = 100;
  int version = 1;
  int fftw_opt = FFTW_ESTIMATE;
  bool simple_parallelism = false;

  parse_arguments(argc, argv, &n, &k, &num_inputs, &version, &fftw_opt,
                  &simple_parallelism);

  std::vector < sfft_plan * >plans;
  if (simple_parallelism)
    {
      for (int i = 0; i < num_inputs; i++)
        {
          sfft_plan *p = plan_by_params(n, k, version, fftw_opt);
          if (p == NULL)
            {
              usage(argv[0]);
              exit(1);
            }
          plans.push_back(p);
        }
    }
  else
    {
      sfft_plan *plan;
      if ((plan = plan_by_params(n, k, version, fftw_opt)) == NULL)
        {
          usage(argv[0]);
          exit(1);
        }
      plans.push_back(plan);
    }

  complex_t **x = generate_input(n, k, num_inputs);
  sfft_output *out = new sfft_output[num_inputs];

  CALLGRIND_START_INSTRUMENTATION;

  timespec ts, te;
  clock_gettime(CLOCK_REALTIME, &ts);
  if (simple_parallelism)
    {
      #pragma omp parallel for
      for (int i = 0; i < num_inputs; ++i)
        {
          sfft_exec(plans[i], x[i], out + i);
        }
    }
  else
    {
      sfft_exec_many(plans[0], num_inputs, x, out);
    }
  clock_gettime(CLOCK_REALTIME, &te);

  CALLGRIND_STOP_INSTRUMENTATION;

  double t =
    (te.tv_sec + 1e-9 * te.tv_nsec) - (ts.tv_sec + 1e-9 * ts.tv_nsec);
  std::cout << "TIME " << t << std::endl;
}
