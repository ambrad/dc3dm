/* dc3dm: Software to form and apply a 3D DDM operator for a nonuniformly
 *        discretized rectangular fault
 *   Version 0.3
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 * https://pangea.stanford.edu/research/CDFM/software
 * dc3dm is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.0
 */
 
This package provides:
- A program to
  - Nonuniformly mesh a rectangular fault in an elastic whole or half space.
  - Calculate the matrix entries of the displacement discontinuity method (DDM)
    linear operator for this nonuniform mesh. This DDM uses Y. Okada's routines;
    see below for the refererence.
  - Create an H-matrix approximation to this matrix.
- An optional Matlab interface to create the input files to build the H-matrix.
- An optional Matlab interface to analyze mesh data.
- An optional Matlab interface to compute MVP.
- C++ routines to compute matrix-vector products (MVP) and related operations
  using this H-matrix approximation.
- A limited C interface to the MVP routines.
- A limited Fortran interface to the MVP routines.
- A limited C++ MPI interface to the MVP routines.

dc3dm is a layer on top of, and comes with, the separate package hmmvp. You
should download and use hmmvp by itself rather than dc3dm if you have your own
Green's function or if your fault model is something more complicated than a
rectangle.


Installation
------------
- If you are working in a Unix-like developer environment, follow these
  instructions. I don't have anything for Windows or Macs that lack Unix-like
  developer tools. But you can still follow the Makefile as a guide to how to
  write the appropriate build script for your environment.
- Edit the top of Makefile with compiler and LAPACK/BLAS information. Set the
  optimization level (opt =) and serial/parallel mode (mode =, with
  (s)erial, OpenMP (omp), and MPI (mpi)).
- On the command line, type 'make'.
- (Optional.) Start Matlab and cd to 'matlab/'. Type 'make'. It seems Macs have
  a problem with recent versions of Matlab; follow these instructions:
      http://www.mathworks.com/matlabcentral/answers/94092


Usage
-----
The basic usage pattern is to write key-value files that describe the problem,
run bin/dc3dm 'mesh', 'build', and 'compress' commands, and then to compute
matrix-vector products with the resulting matrix in your own code by linking
against lib/libhmmvp.

Instructions and examples:
- In a terminal, type
      ./bin/dc3dm help
      ./bin/dc3dm help mesh
      ./bin/dc3dm help build
      ./bin/dc3dm help compress
- In Matlab, type
      addpath matlab; help dc3dm; help hmmvp;
- See examples/ex.m for Matlab usage. ('addpath examples/' to access it in
  Matlab.) Even if you don't intend to use dc3dm in Matlab, reading through ex.m
  can be useful.
- As explained in greater detail in ex.m and the dc3dm help, you need to create
  key-value files to describe your problem. This can be done in Matlab using
  matlab/dc3dm.Read/WriteKvf (no mex file needed) or in C++ using
  util/include/KeyValueFile.hpp
- See examples/mvp_*.c* for example matrix-vector product usage.
- To apply an H-matrix, link against lib/libhmmvp_*.a, lapack, and blas. Follow
  the example in the Makefile target 'mvp'.

The Matlab interface makes using dc3dm simpler. However, it is not
necessary. util/include/KeyValueFile.hpp creates and reads key-value files in
C++. The .elem file emitted by the 'build' operation provides the necessary
information about the mesh and element ordering.


If you use this software in research that you publish, I would appreciate your
citing this software package as follows:

    A. M. Bradley, Software for efficient static dislocation-traction
    calculations in fault simulators, Seis. Res. Lett. 85(6), 2014.

Feel free to contact me at (my permanent email address) ambrad@cs.stanford.edu
to check for a new version of this software. In addition, current versions are
maintained on these sites:
    pangea.stanford.edu/research/CDFM/software
    cs.stanford.edu/~ambrad


Limitations
-----------
- Not all transpose MVP functionality is implemented.
- bin/dc3dm is parallelized by OpenMP but not MPI; however, libhmmvp uses MPI if
  desired (for MVP).


Plans
-----
- Command 'mesh-rsf' that implements one example of a resolution function.
- Command 'compress-jxk' to implement multicomponent problems using a tolerance
  appropriate for the whole problem.
- Command 'est-con-tol' that estimates a conservative tolerance.
- A python interface.


Acknowledgements
----------------
Paul Segall, my postdoc advisor, and I gratefully acknowledge support from the
National Science Foundation (EAR-0838267), the U.S. Geological Survey
(G12AP20030), and SCEC (13096). I would also like to thank Paul for his support
of this work. (Any errors are mine.)

dc3dm uses the code dc3.f (with modifications described at the top of
external/dc3omp.f) of Y. Okada:
    Y. Okada, Internal deformation due to shear and tensile faults in a
    half-space, Bull. Seism. Soc. Am., 1992.
