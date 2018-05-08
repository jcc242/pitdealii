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
That time-dependent data will be provided at each timestep by XBraid via the <kbd>vector</kbd> struct.

### Math details

In the default mode, this code solves the heat equation with,

$$\frac{\partial u}{\partial t} - \Delta u = f(\boldsymbol{x},t), \qquad &\forall\boldsymbol{x}\in\Omega,t\in\left( 0,T \right),$$
with initial conditions,
$$u(\boldsymbol{x},0) = u_0(\boldsymbol{x}), \qquad &\forall \boldsymbol{x}\in\Omega,$$
where $u_0(\boldsymbol{x})=0$, and boundary conditions,
$$u(\boldsymbol{x},t) = g(\boldsymbol{x},t), \qquad &\forall \boldsymbol{x}\in\partial\Omega,t\in\left( 0,T \right)$$
where,
$$g(\boldsymbol{x},t) = 0,$$
and forcing function,
$$
  f(\mathbf x, t)
  =
  \left\{
  \begin{array}{ll}
    \chi_1(\mathbf x) & \text{if  \(x>0.5\) and \(y>-0.5\)}
    \\
    \chi_2(\mathbf x)
    & \text{if \(x>-0.5\) and \(y>0.5\)}
    \\
    0 & \text{otherwise}
  \end{array}
  \right.
$$
with,
$$
	\chi_1(\mathbf x) = \exp{-0.5\frac{(time-0.125)^2}{0.005}}
$$
and,
$$
  \chi_2(\mathbf x) = \exp{-0.5\frac{(time-0.375)^2}{0.005}}
$$
The forcing function is a Gaussian pulse in time that is centered around 0.125 time units for $\chi_1$ and 0.375 time units for $\chi_2$.
A gaussian function was chosen because it is a continuous function in time and so the solution state can be compared bitwise with the serial-in-time solution.

The method of manufactured solutions is used to test the correctness of the implementation.
In the method of manufactured solutions, we create a solution $u_h$ to the heat equation, then compute the boundary conditions, initial conditions, and forcing functions required to generate that solution.
The created solution selected is,
$$u_h = e^{-4\pi^2t}cos(2 \pi x)cos(2 \pi y), \qquad \forall \boldsymbol{x} \in \Omega \cup \partial\Omega$$
with derivatives,
$$\frac{\partial u}{\p artial t} &= -4 \pi^2 e^{-4\pi^2t}cos(2 \pi x)cos(2 \pi y),$$
  $$-\Delta u &= 8 \pi^2 e^{-4\pi^2t}cos(2 \pi x)cos(2 \pi y)$$
  $$\frac{\partial u}{\partial x} &= -2 \pi e^{-4\pi^2t}sin(2\pi x)cos(2\pi y)$$
  $$\frac{\partial u}{\partial x} &= -2 \pi e^{-4\pi^2t}cos(2\pi x)sin(2\pi y)$$
and therefore we specify the terms from the governing equations as,
$$f(\boldsymbol{x},t) = 4 \pi^2 e^{-4\pi^2t}\cos(2 \pi x)\cos(2 \pi y), \qquad&\forall\boldsymbol{x}\in\Omega,t\in\left( 0,T \right),$$
$$u_0(\boldsymbol{x}) = cos(2 \pi x)\cos(2 \pi y), \qquad &\forall \boldsymbol{x}\in\Omega,$$
$$g(\boldsymbol{x},t) = e^{-4\pi^2t}\cos(2 \pi x)\cos(2 \pi y), \qquad &\forall \boldsymbol{x} \in \partial\Omega.$$

The manufactured solution is run on progressively more refined grids and the solution generated by the finite element method is compared to the exact solution $u_h$.
The convergence rate of the error is calculated with
\begin{align*}
  \Delta \epsilon_n = \frac{\ln{\epsilon_{n-1}/\epsilon_{n}}}{\ln{r_n}}
\end{align*}
where $\Delta \epsilon_n$ is the convergence rate of error $\epsilon$ between a mesh $n$ and
coarser mesh $n-1$ that have a refinement ratio of $r_n$.
Shown in Table~\ref{tbl:convergenceRate} is the convergence rate as the mesh is
refined.
The $\Delta t$ is reduced by a factor of 2 for every global refinement of the mesh.

\begin{table}\label{tbl:errorMFG}
  \begin{center}
    \begin{tabular}{|c|c|c|c|c|c|} \hline
      cycle & \# cells & \# dofs & $L^2$-error & $H^1$-error & $L^\infty$-error \\ \hline
      125 & 48 & 65 & 6.036e-03 & 6.970e-02 & 7.557e-03\\ \hline
      250 & 192 & 225 & 1.735e-03 & 3.414e-02 & 2.721e-03 \\ \hline
      500 & 768 & 833 & 4.513e-04 & 1.690e-02 & 7.410e-04 \\ \hline
      1000 & 3072 & 3201 & 1.140e-04 & 8.426e-03 & 1.877e-04 \\ \hline
      2000 & 12288 & 12545 & 2.859e-05 & 4.209e-03 & 4.715e-05 \\ \hline
    \end{tabular}
  \end{center}
  \caption{The error of the simulation as the mesh is globally refined.}
\end{table}

\begin{table}\label{tbl:convergenceRate}
  \begin{center}
    \begin{tabular}{|c|c|c|c|c|c|} \hline
      cycle & \# cells & \# dofs & Slope $L^2$ & Slope $H^1$  & Slope $\textrm{L}^\infty$ \\ \hline
      125 & 48 & 65 & --- & --- & --- \\ \hline
      250 & 192 & 225 & 1.798 & 1.030 & 1.474 \\ \hline
      500 & 768 & 833 & 1.943 & 1.014 & 1.877 \\ \hline
      1000 & 3072 & 3201 & 1.985 & 1.004 & 1.981 \\ \hline
      2000 & 12288 & 12545 & 1.995 & 1.001  & 1.993 \\ \hline
    \end{tabular}
  \end{center}
  \caption{Convergence rate of the error from Table~\ref{tbl:errorMFG}.}
\end{table}



### Code Organization

#### The src directory

The entry point of the code is in <kbd>pitdealii.cc</kbd> and sets up XBraid for a simulation. The XBraid setup involves initializing the <kbd>app</kbd> struct and configuring XBraid for the desired number of timesteps, number of iterations, etc.
The functions implemented for XBraid's use are declared in <kbd>BraidFuncs.hh</kbd> and defined in <kbd>BraidFuncs.cc</kbd>. The <kbd>HeatEquation</kbd> class and all <kbd>deal.ii</kbd> functionality is declared in <kbd>HeatEquation.hh</kbd> and defiend in <kbd>HeatEquationImplem.hh</kbd>. Since <kbd>HeatEquation</kbd> is a class template, its definition file <kbd>HeatEquationImplem.hh</kbd> is included at the bottom of <kbd>HeatEquation.hh</kbd>. Lastly various helper functions and variables such as the current processor id and the output stream are declared in <kbd>Utilities.hh</kbd> and defined in <kbd>Utilities.cc</kbd>.

#### The test directory

This directory contains tests to be built and run with CMake.
These tests verify the correct implementation of the various functions.

#### The doc directory

This directory contains some further documentation on the code.
There is a doxygen config file in the <kbd>doxygen/doc/</kbd> folder, you can generate doxygen documentation with <kbd>doxygen pitdealii.doxygen</kbd>.
In the design folder there is a latex file that has some further documentation on the algorithm implemented in this software.

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

  * <kbd>CMAKE_BUILD_TYPE=Debug/Release</kbd>
  * <kbd>DO_MFG=ON/OFF</kbd>
  * <kbd>USE_MPI=ON/OFF</kbd>

The build type specifies whether to compile with debugging symbols, assertions, and optimizations or not.

The option for manufactured solutions (<kbd>DO_MFG</kbd>) switches from solving the "standard" heat equation to solving a known solution heat equation so that the correctness of the code can be tested.

Lastly the MPI option is only used to specify where to write output information when using the <kbd>pout()</kbd> function from the Utilities.hh file.
If <kbd>USE_MPI</kbd> is set to <kbd>ON</kbd>, then every processor writes to its own file called <kbd>pout.<#></kbd> where <kbd><#></kbd> is the processor number.
If <kbd>USE_MPI</kbd> is set to <kbd>OFF</kbd>, then every processor writes to stdout.

## Running

Once pitdealii has been compiled, the program can be run by calling the binary generated in <kbd>./build/bin/</kbd>.
The test can be run by calling <kbd>ctest</kbd> from inside the <kbd>build/</kbd> directory.
Unless the output path has been changed in the source code (currently hardcoded), then the output files will be placed into the folder the command was called from.

## Results
