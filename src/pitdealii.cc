/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2013 - 2017 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2013
 */

#include <fstream>
#include <iostream>

#include "BraidFuncs.hh"
#include "HeatEquation.hh"

int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;

      /* Initialize MPI */
      MPI_Comm      comm, comm_x, comm_t;
      int rank;
      MPI_Init(&argc, &argv);
      comm   = MPI_COMM_WORLD;
      MPI_Comm_rank(comm, &rank);

      // Set up X-Braid
      /* Initialize Braid */
      braid_Core core;
      double tstart = 0.0;
      double tstop = 0.5;
      double ntime = 500;
      my_App *app = new(my_App);
      
      app->eq.define();

      // HeatEquation<2> heat_equation_solver;

      braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
                 my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
                 my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

      /* Define XBraid parameters 
       * See -help message for descriptions */
      int       max_levels    = 2;
      int       nrelax        = 1;
      int       skip          = 0;
      double    tol           = 1.0e-07;
      int       cfactor       = 2;
      int       max_iter      = 30;
      int       min_coarse    = 3;
      int       fmg           = 0;
      int       scoarsen      = 0;
      int       res           = 0;
      int       wrapper_tests = 0;
      int       print_level   = 1;
      int       access_level  = 1;
      int       use_sequential= 0;

      braid_SetPrintLevel( core, print_level);
      braid_SetAccessLevel( core, access_level);
      braid_SetMaxLevels(core, max_levels);
      braid_SetMinCoarse( core, min_coarse );
      braid_SetSkip(core, skip);
      braid_SetNRelax(core, -1, nrelax);
      braid_SetAbsTol(core, tol);
      braid_SetCFactor(core, -1, cfactor);
      braid_SetMaxIter(core, max_iter);
      braid_SetSeqSoln(core, use_sequential);

      braid_Drive(core);

      braid_Destroy(core);
      // heat_equation_solver.run(tstart, tstop, ntime);

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what()
                << std::endl << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!"
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
