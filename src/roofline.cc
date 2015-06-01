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
#include "measuring_core.h"

int main(int argc, char *argv[])
{
  long counters[] = { 0x10, 0x80, 0x10, 0x10, 0x11, 0x02, 0x11, 0x01 };

  if (0 != measurement_init(counters, 0, 0))
    {
      //      std::cerr << "Failed to initialize measuring core.\n";
      //      return 1;
    }

  simulation sim(argc, argv);
  measurement_start();
  for (int i = 0; i < 10; i++)
    sim.run();
  measurement_stop(1);

  measurement_end();

  return 0;
}
