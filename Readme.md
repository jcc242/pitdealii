# Parallel in Time deal.ii

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
