# Development of a Julia Code for Turbulent Flow Simulations

[Manuel F. Schmid](mailto:mfs2173@columbia.edu) • Marco G. Giometto • Gregory A. Lawrence • Marc B. Parlange

## Motivation & Objectives

When developing a fluid dynamics code for turbulence research, there are three
important goals:

1) **Correctness**: While all simulations aim to be implement the numerical
   methods correctly, this is particularly important in turbulence research.
   Due to the chaotic nature of turbulent flows, simulations are difficult to
   validate. Furthermore, the models are regularly adapted to the problem at
   hand and often novel simulation approaches are evaluated, where it is
   particularly important that the implemented model corresponds exactly to its
   mathematical definition.

2) **Performance**: Turbulence simulations are often limited by their
   computational performance. When single runs on a powerful parallel system
   regularly take weeks, there is little room to trade off performance for
   convenience. Often there is a direct link between the performance of a code
   and its usefulness for research.

3) **Adaptability**: Many applications of computational fluid dynamics require
   problem-specific adaptations of the simulation code, for example to account
   for processes at length scales that are too large or too small to resolve or
   to evaluate a new parameterization.

The pseudo-spectral Fortran code originally developed by John Albertson has
been in use in our academic family for over two decades and has been the
foundation for a lot of useful research. While it has been extensively
validated and its performance has remained satisfactory, making substantial
changes is tedious and error-prone.

The Julia programming language has been developed specifically for scientific
computations. It promises to provide the convenience of a high-level language
similar to Python, MATLAB, and R without sacrificing performance. The language
is open source and has reached version 1.0 in August 2018.

This project attempts to develop a Julia code for turbulent flow simulations
that has the same numerical set-up as the Parlange–Albertson code and achieves
similar performance, but provides more flexibility for future development and
stronger guarantees for its current and future correctness. While the new code
won’t support all functionality of the Parlange–Albertson code initially, it
should provide the dynamical core for direct numerical simulation and basic
large-eddy simulation of turbulent channel flows.

## Data & Code

The repository contains the code and data required to reproduce the manuscript “BoundaryLayerDynamics.jl v1.0: a modern codebase for atmospheric boundary-layer simulations”.
For complete reproduction including the simulation runs, the two submodules with the simulation code need to be initiated (Parlange–Albertson code requires access to private repository).

- To build the manuscript, run `make manuscript`. This relies on the figures, which are included in the repository.
- To recreate the figures, run `make figures`. This relies on the precomputed profile data in the `data` directory (HDF5 files).
- To recompute the profiles, run `make profiles`. This relies on simulation snapshots and performance data in the `data/*` subdirectories. These are not included in the repository due to their larger size, but can be recreated with the available code.
- To rerun validation simulations, run `make dns-julia`, `make les-julia`, or `make les-fortran`. This requires 64 MPI processes (can be reconfigured) and may take several hours.
- To rerun performance tests, run `make perf`. This is meant to be run on the Stampede2 system and has to be adapted when running on a different system or using a different account.

The data was generated with Julia 1.9 but any recent Julia version (≥1.6) should work.
The exact versions of dependencies are given in the `Manifest.toml` files and should be downloaded automatically when using the above `make` commands.

## Goals & Milestones

- [x] Feasibility study: Check whether Julia has all the required functionality
  to run a simple test with FFTs and MPI communication in an HPC setting.
- [x] Basic functionality: Implement the core functionality for DNS of a simple
  channel flow, including MPI parallelization.
- [x] Evaluate performance and attempt to bring it in line with the
  Parlange–Albertson code.
- [x] Validation against laminar flow simulations, with automated testing.
- [x] Full support for DNS with the required functionality for initialization,
  restart, output of snapshots, and flow statistics, validated with the Lee
  & Moser (2015) data.
- [x] LES with a static Smagorinsky SGS model with support for computation of
  spectra, validated against Parlange–Albertson code.
- [x] Publication (internal or public) of code with documentation of the
  numerical approach and guidelines for further development.
- [-] Publication of manuscript describing the methods, validation, and
  performance of the code.

## Results

- The code produces correct solutions to laminar flow problems with the expected convergence properties.
- The code produces DNS results that match those published by Lee & Moser (2015).
- The code produces LES results with the expected profiles & spectra.
- The code matches the performance of the Parlange–Albertson code.

## Publications

- Poster AGU 2018
- Poster Burgers Summer School 2019
- Presentation at APS DFD 2019
- Paper GMD (in progress)
