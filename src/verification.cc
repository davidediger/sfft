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

#include "simulation.h"

int main(int argc, char **argv)
{
  const double NOISE_THRESHOLD = 0.1;
  const double ERROR_THRESHOLD = 0.1;

  simulation sim(argc, argv);

  int n = sim.get_n();
  complex_t *f = sim.frequencies();
  sfft_output ans = sim.run();

  for (int i = 0; i < n; i++)
    {
      if (cabs(f[i]) > NOISE_THRESHOLD)
        {
          // Was the frequency recovered?
          if (ans.find(i) == ans.end())
            {
              std::cout << "ERROR: Frequency " << i
                        << " was not recovered!\n";
              return 1;
            }
          // Was the frequency recovered with the right value?
          if (cabs(ans[i] - f[i]) > ERROR_THRESHOLD)
            {
              std::cout << "ERROR: Error of frequency " << i
                        << "is too big:\n"
                        << "  Expected: " << f[i] << "\n"
                        << "  Actual  : " << ans[i] << "\n";
              return 2;
            }
        }
    }

  std::cout << "OK\n";
  return 0;
}
