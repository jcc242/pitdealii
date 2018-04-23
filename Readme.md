# Parallel in Time deal.ii

## Overview

Over the last few years, the clock speed per processor core has stagnated.
This stagnation has lead to the design of larger core counts in high performance computing machines.
As a result of these developments, increased concurrency in numerical algorithms must be developed in order to take advantage of this architecture style.
Perhaps the largest bottleneck in concurrency for time-dependent simulations is traditional time integration methods.
Traditional time integration methods solve for the time domain sequentially, and as the spatial grid is refined a proportionally larger number of time steps must be taken to maintain accuracy and stability constraints.
While solving the time domain sequentially with a traditional time integration method is an optimal algorithm of order $\mathcal{O}(n)$, the $n$ time steps are not solved concurrently.

The goal of this project is to make use of the XBraid library from Lawrence Livermore National Laboratory to solve the time domain in parallel using multigrid-reduction-in-time techniques.
The XBraid library is implemented in C and aims to be a non-intrusive method to implement parallel time marching methods into existing codes.

## Implementation

### XBraid introduction

In order to use the XBraid library, several data structures and functions must be implemented and provided to the XBraid solver struct.
The two required data structures are the <kbd>app</kbd> and <kbd>vector</kbd> structures.
In general, the <kbd>app</kbd> struct contains the time independent data and the <kbd>vector</kbd> struct contains the time dependent data.
For this initial example, the time independent data includes the mesh which is fixed for all time steps, and the time dependent data is the solution state vector.
The functions tell XBraid how to perform operations on the data type used by your solver, in this case deal.ii uses the Vector<double> data type.
These operations include how to initialize the data at a given time, how to sum the data, and how to pack and unpack linear buffers for transmission to other processors via MPI.
The XBraid documentation should be read for a full list of functions that must be implemented and the details of what the function should do.
The typical format is the function is called with arguments of the <kbd>app</kbd> struct, one or more <kbd>vector</kbd> structs, and a <kbd>status</kbd> struct that contains information on the current status of the XBraid simulation (the current multigrid iteration, the level the function is being called from, the time and timestep number, etc.).

Perhaps the most important function is the <kbd>step</kbd> function.
This function tells XBraid how to advance the solution forward in time from the initial to the final times given in the <kbd>status</kbd> struct.
This method uses a traditional time integration method such as the fourth order explicit Runge Kutta method.

### deal.ii details

The solver used in this example is based off the heat equation solver from the [Step-26 Tutorial](https://dealii.org/developer/doxygen/deal.II/step_26.html "Step-26 Tutorial").
The <kbd>HeatEquation</kbd> class becomes member data to XBraid's <kbd>app</kbd> struct, and XBraid's <kbd>vector</kbd> struct becomes a wrapper for deal.ii's <kbd>Vector<double></kbd> data type.
The <kbd>HeatEquation</kbd> class cannot simply be used as is though as it contains both time dependent and time independent member data.
In order to simplify the problem the adaptive mesh refinement is removed.
Theoretically XBraid is capable of working with adaptive mesh refinement and in fact contains support for time refinement (which is also not used for simplicity).
All adaptive mesh refinement functionality is removed from the sovler.
The time-dependent solution state vectors are also removed from the <kbd>HeatEquation</kbd> member data.
That data will be provided at each timestep by XBraid via the <kbd>vector</kbd> struct.

## Compiling

To compile, you need deal.ii and XBraid to be installed with development headers somwehere on your system.
Some implementation of MPI such as OpenMPI with development headers must also be installed.
The source code for deal.ii is available at [https://dealii.org/](https://dealii.org/) and the source code for XBraid is available at [https://computation.llnl.gov/projects/parallel-time-integration-multigrid](https://computation.llnl.gov/projects/parallel-time-integration-multigrid).
See the documentation of each package for compilation and installation instructions.

Depending on where they are installed, pitdealii may need help finding these libraries.
To find deal.ii, pitdealii first looks in typical deal.ii install directories followed by one directory up(<kbd>../</kbd>), two directories up (<kbd>../../</kbd>), and lastly in the environment variable <kbd>DEAL_II_DIR</kbd>.
In contrast, XBraid currently does not have any default locations to look for and so the environment variable <kbd>BRAID_DIR</kbd> must be specified.
For MPI, pitdealii looks in standard installation folders only, for that reason I recommend you install MPI with your package manager.

A compile process of pitdealii may look like,

    mkdir build
	cd build
	BRAID_DIR=/path/to/braid/ cmake ../
	make

There is currently no option to install pitdealii anywhere.
The binaries are generated in the <kbd>bin</kbd> folder, and tests are placed into the <kbd>test</kbd> folder.
Options that can be passed to CMake for pitdealii include:

  * CMAKE_BUILD_TYPE=Debug/Release
  * DO_MFG=ON/OFF
  * USE_MPI=ON/OFF

The build type specifies whether to compile with debugging symbols, assertions, and optimizations or not.

The option for manufactured solutions (<kbd>DO_MFG</kbd>) switches from solving the "standard" heat equation to solving a known solution heat equation so that the correctness of the code can be tested.

Lastly the MPI option is only used to specify where to write output information when using the <kbd>pout()</kbd> function from the Utilities.hh file.
If <kbd>USE_MPI</kbd> is set to <kbd>ON</kbd>, then every processor writes to its own file called <kbd>pout.<#></kbd> where <kbd><#></kbd> is the processor number.
If <kbd>USE_MPI</kbd> is set to <kbd>OFF</kbd>, then every processor writes to stdout.

## Running

Once pitdealii has been compiled, the program can be run by calling the binary generated in <kbd>./build/bin/</kbd>.
The test can be run by calling <kbd>ctest</kbd> from inside the <kbd>build/</kbd> directory.
Unless the output path has been changed in the source code (currently hardcoded), then the output files will be placed into the folder the command was called from.
