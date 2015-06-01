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

#include "parameters.h"

void
get_expermient_vs_N_parameters(int N, bool WITH_COMB, double &Bcst_loc,
                               double &Bcst_est, double &Comb_cst,
                               int &loc_loops, int &est_loops,
                               int &threshold_loops, int &comb_loops,
                               double &tolerance_loc, double &tolerance_est)
{

  if (WITH_COMB)
    {

      switch (N)
        {
        case 8192:
          Bcst_loc = 2;
          Bcst_est = 2;
          Comb_cst = 32;
          comb_loops = 8;
          est_loops = 16;
          loc_loops = 7;
          threshold_loops = 6;
          tolerance_loc = 1e-8;
          tolerance_est = 1e-8;
          break;
        case 16384:
          Bcst_loc = 4;
          Bcst_est = 4;
          Comb_cst = 32;
          comb_loops = 8;
          est_loops = 10;
          loc_loops = 6;
          threshold_loops = 5;
          tolerance_loc = 1e-8;
          tolerance_est = 1e-8;
          break;
        case 32768:
          Bcst_loc = 4;
          Bcst_est = 2;
          Comb_cst = 64;
          comb_loops = 4;
          est_loops = 8;
          loc_loops = 5;
          threshold_loops = 4;
          tolerance_loc = 1e-8;
          tolerance_est = 1e-8;
          break;
        case 65536:
          Bcst_loc = 4;
          Bcst_est = 2;
          Comb_cst = 128;
          comb_loops = 6;
          est_loops = 10;
          loc_loops = 4;
          threshold_loops = 2;
          tolerance_loc = 1e-8;
          tolerance_est = 1e-8;
          break;
        case 131072:
          Bcst_loc = 1;
          Bcst_est = 1;
          Comb_cst = 8;
          comb_loops = 2;
          est_loops = 12;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 262144:
          Bcst_loc = 1;
          Bcst_est = 1;
          Comb_cst = 8;
          comb_loops = 2;
          est_loops = 14;
          loc_loops = 5;
          threshold_loops = 4;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 524288:
          Bcst_loc = 0.5;
          Bcst_est = 0.5;
          Comb_cst = 8;
          comb_loops = 1;
          est_loops = 10;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 1048576:
          Bcst_loc = 0.5;
          Bcst_est = 0.5;
          Comb_cst = 8;
          comb_loops = 2;
          est_loops = 12;
          loc_loops = 4;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 2097152:
          Bcst_loc = 0.5;
          Bcst_est = 0.2;
          Comb_cst = 8;
          comb_loops = 1;
          est_loops = 10;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 4194304:
          Bcst_loc = 0.5;
          Bcst_est = 0.2;
          Comb_cst = 8;
          comb_loops = 1;
          est_loops = 8;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 8388608:
          Bcst_loc = 0.5;
          Bcst_est = 0.2;
          Comb_cst = 8;
          comb_loops = 1;
          est_loops = 8;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 16777216:
          Bcst_loc = 0.5;
          Bcst_est = 0.2;
          Comb_cst = 16;
          comb_loops = 1;
          est_loops = 8;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        }
    }
  else
    {
      switch (N)
        {
        case 8192:
          Bcst_loc = 2;
          Bcst_est = 2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 16;
          loc_loops = 7;
          threshold_loops = 6;
          tolerance_loc = 1e-8;
          tolerance_est = 1e-8;
          break;
        case 16384:
          Bcst_loc = 4;
          Bcst_est = 4;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 10;
          loc_loops = 6;
          threshold_loops = 5;
          tolerance_loc = 1e-8;
          tolerance_est = 1e-8;
          break;
        case 32768:
          Bcst_loc = 4;
          Bcst_est = 2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 8;
          loc_loops = 5;
          threshold_loops = 4;
          tolerance_loc = 1e-8;
          tolerance_est = 1e-8;
          break;
        case 65536:
          Bcst_loc = 4;
          Bcst_est = 2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 8;
          loc_loops = 5;
          threshold_loops = 4;
          tolerance_loc = 1e-8;
          tolerance_est = 1e-8;
          break;
        case 131072:
          Bcst_loc = 2;
          Bcst_est = 1;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 10;
          loc_loops = 5;
          threshold_loops = 4;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 262144:
          Bcst_loc = 2;
          Bcst_est = 0.5;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 14;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 524288:
          Bcst_loc = 1;
          Bcst_est = 0.5;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 12;
          loc_loops = 5;
          threshold_loops = 4;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 1048576:
          Bcst_loc = 2;
          Bcst_est = 0.5;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 12;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 2097152:
          Bcst_loc = 2;
          Bcst_est = 0.2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 15;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 4194304:
          Bcst_loc = 4;
          Bcst_est = 0.2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 10;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 8388608:
          Bcst_loc = 2;
          Bcst_est = 0.2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 8;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;
        case 16777216:
          Bcst_loc = 4;
          Bcst_est = 0.2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 8;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1e-8;
          break;

        }
    }

  return;
}

void
get_expermient_vs_K_parameters(int K, bool WITH_COMB, double &Bcst_loc,
                               double &Bcst_est, double &Comb_cst,
                               int &loc_loops, int &est_loops,
                               int &threshold_loops, int &comb_loops,
                               double &tolerance_loc, double &tolerance_est)
{

  if (WITH_COMB)
    {

      switch (K)
        {
        case 50:
          Bcst_loc = 0.5;
          Bcst_est = 0.2;
          Comb_cst = 16;
          comb_loops = 1;
          est_loops = 10;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1.0e-8;
          break;
        case 100:
          Bcst_loc = 0.5;
          Bcst_est = 0.2;
          Comb_cst = 16;
          comb_loops = 1;
          est_loops = 12;
          loc_loops = 4;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1.0e-8;
          break;
        case 200:
          Bcst_loc = 0.5;
          Bcst_est = 0.5;
          Comb_cst = 32;
          comb_loops = 1;
          est_loops = 8;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-6;
          tolerance_est = 0.5e-8;
          break;
        case 500:
          Bcst_loc = 0.5;
          Bcst_est = 0.5;
          Comb_cst = 64;
          comb_loops = 1;
          est_loops = 10;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-6;
          tolerance_est = 0.5e-8;
          break;
        case 1000:
          Bcst_loc = 1;
          Bcst_est = 1;
          Comb_cst = 128;
          comb_loops = 3;
          est_loops = 12;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-6;
          tolerance_est = 0.5e-8;
          break;
        case 2000:
          Bcst_loc = 1;
          Bcst_est = 1;
          Comb_cst = 512;
          comb_loops = 3;
          est_loops = 16;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-7;
          tolerance_est = 0.2e-8;
          break;
        case 2500:
          Bcst_loc = 1;
          Bcst_est = 1;
          Comb_cst = 512;
          comb_loops = 3;
          est_loops = 16;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-7;
          tolerance_est = 0.2e-8;
          break;
        case 4000:
          Bcst_loc = 1;
          Bcst_est = 2;
          Comb_cst = 512;
          comb_loops = 3;
          est_loops = 14;
          loc_loops = 8;
          threshold_loops = 7;
          tolerance_loc = 1e-8;
          tolerance_est = 0.5e-8;
          break;
        }
    }
  else
    {
      switch (K)
        {
        case 50:
          Bcst_loc = 4;
          Bcst_est = 0.2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 10;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1.0e-8;
          break;
        case 100:
          Bcst_loc = 2;
          Bcst_est = 0.2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 12;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 1.0e-8;
          break;
        case 200:
          Bcst_loc = 4;
          Bcst_est = 0.5;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 10;
          loc_loops = 3;
          threshold_loops = 2;
          tolerance_loc = 1e-6;
          tolerance_est = 0.5e-8;
          break;
        case 500:
          Bcst_loc = 2;
          Bcst_est = 1;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 12;
          loc_loops = 4;
          threshold_loops = 3;
          tolerance_loc = 1e-6;
          tolerance_est = 0.5e-8;
          break;
        case 1000:
          Bcst_loc = 2;
          Bcst_est = 1;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 12;
          loc_loops = 5;
          threshold_loops = 4;
          tolerance_loc = 1e-6;
          tolerance_est = 1.0e-8;
          break;
        case 2000:
          Bcst_loc = 2;
          Bcst_est = 1;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 16;
          loc_loops = 5;
          threshold_loops = 4;
          tolerance_loc = 1e-7;
          tolerance_est = 0.5e-8;
          break;
        case 2500:
          Bcst_loc = 2;
          Bcst_est = 1;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 16;
          loc_loops = 5;
          threshold_loops = 4;
          tolerance_loc = 1e-7;
          tolerance_est = 0.5e-8;
          break;
        case 4000:
          Bcst_loc = 2;
          Bcst_est = 2;
          Comb_cst = 1;
          comb_loops = 1;
          est_loops = 14;
          loc_loops = 6;
          threshold_loops = 5;
          tolerance_loc = 1e-8;
          tolerance_est = 1.0e-8;
          break;
        }
    }

  return;
}
