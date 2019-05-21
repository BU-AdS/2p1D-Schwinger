# (2+1)D-Schwinger Model

C++ code for a 2D Schwinger Model, with the option to add a 3rd dimension
for the gauge field with open boundary conditions and unit valued links in
the extra dimension. Fermions are restricted to the central slice.

## Features

The code is 'built for comfort, not for speed'. No attempts have been made to
optimise the code for performance. Each fermion action (and each instance in
2D or 3D) is given its own directory with its own main.cpp file and Makefile.

### Actions

We use the Wilson gauge action throughout. At the present time, we offer

   1. 2 flavour Staggered
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
   5. Vacuum trace
   5. Pion correlation function

### Usage

In each (Action/Dimension) directory we provide four files:

   Makefile_template
   main_template.cpp
   looper.sh
   launcher.sh.

We have deliberatly left some parameters of the code as preprocessor defines.
This is so that any attempt to parallelise the code with, say OpenACC, will
be easier. In order to constrcut an exectuable, one must edit the `launcher.sh`
file with the desired parameters. The `launcher.sh` script will then construct
a `Makefile`, a `main.cpp` file, and an executable, and will then launch the job.

## Dependencies

The sole dependency is from ARPACK and is entirely optional. We have tested
against macOS Mojave and used `brew install arpack` to build ARPACK, and link
accordingly.
