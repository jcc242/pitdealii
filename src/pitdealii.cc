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
#include "Utilities.hh"

int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;

      /* Initialize MPI */
      MPI_Comm      comm; //, comm_x, comm_t;
      int rank;
      MPI_Init(&argc, &argv);
      comm   = MPI_COMM_WORLD;
      MPI_Comm_rank(comm, &rank);
      procID = rank;

      // Set up X-Braid
      /* Initialize Braid */
      braid_Core core;
      double tstart = 0.0;
      double tstop = 0.5;
      double ntime = 2500;
      my_App *app = new(my_App);


      // Set up the time and space MPI comms
      // braid_SplitCommworld(&comm, 1, &comm_x, &comm_t);


      // HeatEquation<2> heat_equation_solver;

      braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
                 my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
                 my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

      /* Define XBraid parameters
       * See -help message forf descriptions */
      int       max_levels    = 2;
      // int       nrelax        = 1;
      // int       skip          = 0;
      double    tol           = 1.e-7;
      // int       cfactor       = 2;
      int       max_iter      = 5;
      // int       min_coarse    = 10;
      // int       fmg           = 0;
      // int       scoarsen      = 0;
      // int       res           = 0;
      // int       wrapper_tests = 0;
      int       print_level   = 1;
      int       access_level  = 1;
      // int       use_sequential= 1;

      braid_SetPrintLevel( core, print_level);
      braid_SetAccessLevel( core, access_level);
      braid_SetMaxLevels(core, max_levels);
      //       braid_SetMinCoarse( core, min_coarse );
      //       braid_SetSkip(core, skip);
      //       braid_SetNRelax(core, -1, nrelax);
      braid_SetAbsTol(core, tol);
      //       braid_SetCFactor(core, -1, cfactor);
      braid_SetMaxIter(core, max_iter);
      // braid_SetSeqSoln(core, use_sequential);

      app->eq.define();
      app->final_step = ntime;

      braid_Drive(core);

      // Free the memory now that we are done
      braid_Destroy(core);

      delete app;

      // Clean up MPI
      // MPI_Comm_free(&comm);
      MPI_Finalize();
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

