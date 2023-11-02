# Overview of Data

Large data files are not included in this repository, but each directory
contains a Makefile that runs the necessary commands to generate the data.
Below is an overview of the data along with instructions to rerun all
simulation/processing steps.


## Performance Measurements

The data in the `performance` directory has been collected on the TACC
Stampede2 cluster.

To repeat the tests, run `make launch-skx` for the Skylake nodes or `make
launch-icx` for the Ice Lake nodes, or simply `make` for both. Individual sizes
can be launched with e.g. `make launch-icx-1280`.

Note: The tests of the Fortran code can crash if two simulations running
simultaneously are trying to access the output files at the same time. Running
them one after the other is recommended, especially for the cases with high
core counts, which run faster and therefore access the output files with
a higher frequency.

Measurements were taken with commit `e8b5562` of the Fortran code, and
`d4de5b8` of the Julia code.
