# (2+1)D-Schwinger Model			23rd February 2019

C++ code for a 2D Schwinger Model, with the option to add a 3rd dimension
for the gauge field with open boundary conditions. Fermions are restricted
to the central slice.

## Features and Usage

The code is 'built for comfort, not for speed'. No attempts have been made to
optimise the code for performance, nor have has any reasonable effort been
made to consolidate repetitious functions or code blocks. Each fermion action
(and each instance in 2D or 3D) is given its own directory with its own
linear algebra helper functions, main.cpp file, external library interface,
and Makefile. Code optimisations are a low priority, but code consolidation
will occur when development permits.

### Actions

We use the Wilson gauge action throughout. At the present time, we offer

   1. 2 flavour staggered
   2. 2 flavour Wilson

fermion actions.

### Simulation

We use Leapfrog HMC throughout the code. You have the option to perform dynamic
or quenched simulations, dictated from the command line.

### Utilities

Each action has several utilities at your disposal. Please refer to the comments in
the source for further information. We list here the features as a synopsis:

    1. Gauge field saving/loading.
    2. Gauge field smearing (APE)
    3. Linear algebra suite (linAlgHelpers.h)
    4. ARPACK driven eigensolver

### Measurements

At the preset time, only gauge measurements are implemented.

   1. Average plaquette
   2. Wilson loops (for Creutz ratios)
   3. Polyakov loops
   4. Topological charge

### Usage

In each action directory we provide a `test.sh` file that gives an exhaustive
list of inputs with description. Adjust the variables in the `test.sh` file
then execute `./test.sh` to run the simulation with the given parameters.

## Building and Dependencies

The sole dependency is from ARPACK. We have tested against macOS Mojave and used
`brew install arpack` to build ARPACK.

When making a `main.cpp` file, specifiy ARPACK linking in the `Makefile` and adjust
the path accordingly. Then type

    `make`

to build.

