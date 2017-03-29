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
 
#ifndef INCLUDE_DC3DM_HELP
#define INCLUDE_DC3DM_HELP

// Some routines and definitions for use in the dc3dm main programs.

/* Help text for programs. I write the text in a plain text file and then run
   this routine:

    function WrapPrintHelp (fn)
      fid = fopen(fn, 'r');
      state = 0;
      while (true)
        ln = fgetl(fid);
        if (~ischar(ln)) break; end
        if (numel(ln) >= 10 && strcmp(ln(1:7), '#define'))
          if (state == 1) fprintf('"\n\n'); end
          fprintf('%s \\\n"', ln);
          state = 1;
        else fprintf('%s\\n\\\n', ln); end
      end
      fprintf('"\n\n');
      fclose(fid);
    end
*/

#define _dc3dm_Header_text_ \
"dc3dm version 0.0beta: Software to\n\
  . nonuniformly (or uniformly) discretize a rectangular planar fault in a\n\
    halfspace,\n\
  . compute an H-matrix approximation to the matrix of slip-stress elasticity\n\
    Green's functions computed using Y. Okada's DC3D,\n\
  . and compute matrix-vector products with this matrix.\n\
A.M. Bradley ambrad@cs.stanford.edu\n\
CDFM Group, Geophysics, Stanford University\n\
\n\
\n\
"

#define _dc3dm_all_PrintHelp_text_ \
"dc3dm\n\
  This message.\n\
\n\
dc3dm <key-value file>\n\
    The key-value file describes the problem. The easiest way to create a\n\
  key-value file is to use matlab/dc3dm.Read/WriteKvf. You can create one in C++\n\
  by using the interface described in util/include/KeyValueFile.hpp and linking\n\
  against libhmmvp.\n\
    The key-value file must contain at least this field:\n\
  command [string]: One of\n\
      mesh: Mesh the rectangle according to a resolution function.\n\
      build: Build the fault.\n\
      compress: Form and compress an elasticity matrix.\n\
    Type 'dc3dm help <command>' for help on each command.\n\
\n\
dc3dm help <command>\n\
\n\
\n\
"

#define _dc3dm_mesh_PrintHelp_text_ \
"mesh\n\
\n\
  mesh_write_filename [string]: The suffix .rmu is appended to it. This file\n\
    contains data about the mesh. Let Nr be the number of rectangles. Call the\n\
    ordering of mesh data 'mesh order'.\n\
\n\
  min_len, max_len [scalar]: Minimum and maximum side length of a\n\
    rectangle. max_len is strictly enforced, but some elements may be smaller\n\
    than min_len, though not not smaller than min_len/2.\n\
\n\
  x, y [1D array]: Tensor mesh points relative to an arbitrary origin. The\n\
    domain extends from x(1) to x(end), y(1) to y(end), and so has the\n\
    dimensions\n\
        xlen = x(end) - x(1)\n\
        ylen = y(end) - y(1).\n\
    We use the word 'square' to mean the rectangle that has side lengths\n\
        xlen / ceil(xlen/max_len)  and  ylen / ceil(ylen/max_len).\n\
    Hence to get geometric squares, xlen and ylen should be divisible by\n\
    max_len. y here is called eta in dc3dmBuild, where a dip angle, and so 3D\n\
    coordinates, become relevant.\n\
\n\
  f [2D array]: ny x nx array of f(x, y). f should be sampled enough that linear\n\
    interpolation is sufficient to represent the continuous function. Then the\n\
    mesh is refined so that max f in a square is equal to the square's side\n\
    length\n\
\n\
  do_uniform [scalar, optional, default 0]: Make the mesh uniform. If do_uniform\n\
    = 1, f is ignored and max_len determines the meshing.\n\
\n\
\n\
"

#define _dc3dm_build_PrintHelp_text_ \
"build\n\
\n\
  mesh_read_filename [string]: .rmu file, without suffix, from dc3dmMesh.\n\
\n\
  build_write_filename [string]: Six files are created with suffixes .rmu, .ra,\n\
    .bmb, .hd, .build, and .elem appended.\n\
      The .rmu file is a revised version of the .rmu file from the 'mesh'\n\
    operation.\n\
      The .rmu, .ra, .bmb, .hd, and .build files are for internal use.\n\
      The .elem file is a plain text file in comma-separated value format. It\n\
    describes the elements. Parse this file to determine how to order a vector\n\
    for the matrix-vector product. Column 1 is an index from 1 to Nr, where Nr\n\
    is the number of rectangular elements. Columns 2 and 3 list the (x, y)\n\
    coordinate of an element's center. Columns 4 and 5 list the along-x and\n\
    along-y lengths of the element.\n\
\n\
  do_fullspace [scalar, optional, default false]: Do the fullspace problem\n\
    rather than the halfspace one if do_fullspace = 1.\n\
\n\
  depth_min [scalar, required if halfspace problem]: Minimum depth of fault.\n\
\n\
  dipdeg [scalar, required if halfspace problem]: Dipping angle of the fault in\n\
    degrees. The fault dips in the N-S direction. If positive, the S side is\n\
    deeper than the N side, and opposite if negative.\n\
\n\
  neighborhood [scalar]: Specify the size of the neighborhood for calculating\n\
    the high-order BEM elements. 0 gives constant-slip elements, 1 makes all\n\
    elements have a neighborhood including at least adjacent elements, and\n\
    larger numbers improve the approximation. hmmvpCompress takes longer with\n\
    larger neighborhood. I recommend 4, which is the default if left\n\
    unspecified.\n\
\n\
  Boundary conditions (BC). There are four types of BC:\n\
      Velocity (v-BC). A large slab adjacent to the boundary imposes a constant\n\
    slip rate.\n\
      0-velocity (0v-BC). The slab is at 0 slip rate. This is treated specially\n\
    for computational reasons.\n\
      Periodic. Opposite sides are mapped to each other.\n\
      Free. This occurs on the boundary that is at 0 depth (depth_min = 0), if\n\
    there is one.\n\
      BC are handled by specifying a set of strings. These are as\n\
    follows. First, the cardinal directions of the fault are\n\
        [e]ast: positive-x (strike) direction\n\
        [n]orth: positive-eta (dip) direction.\n\
    [w]est and [s]outh follow from these. Specify a subset of these keys:\n\
        [e|n|s|w]0vbc [scalar]: 0v-BC on one of the cardinal sides.\n\
        [e|n|s|w]vbc [scalar]: v-BC on one of the cardinal sides.\n\
        [ew|ns]pbc [scalar]: Priodic on the east-west or north-south sides.\n\
    The default for all sides is [e|n|w|s]0vbc set to 0.\n\
      The value associated with one of these keys is the dominance of the\n\
    boundary condition; a higher number increases dominance. For example, if\n\
    this is specified:\n\
        evbc: 1, wvbc: 2, n0vbc: 3,\n\
    then (1) 's0vbc' is assumed to have been specified with value 0 and (2)\n\
    large slabs used to approximate the v-BC extend deep south (east and west\n\
    are dominant relative to south) but end at the north side (since the north\n\
    boundary is specified as dominant). Here is a diagram of this example:\n\
                               N: 0 v-BC\n\
                  __________________________________\n\
                           |              | \n\
                   W: v-BC |    domain    | E: v-BC\n\
                           |______________|\n\
                           |   S: 0 v-BC  |\n\
    The issue of ordering matters only for (0-)velocity BC; values for periodic\n\
    BC are ignored. If two values are the same and the boundaries interfere, the\n\
    results are undefined.\n\
      If depth_min is 0, then the boundary at the surface has a free boundary\n\
    condition. This overrides any specification.\n\
      Generally, E-W (along-stride) periodic BC make sense but N-S ones do not\n\
    because of the halfspace.\n\
      As a hint, a subducting fault going to the surface, with y = 0 on the\n\
    updip end, is specified with dipdeg < 0, depth_min = 0. disl_dip > 0 is a\n\
    dislocation in the +y direction. Quite likely, a velocity BC is specified on\n\
    the N side: nvbc = 1; a free surface is on the S side, so it does not have\n\
    to be specified; and the E-W sides are one of these: (1) 0-velocity BC,\n\
    which hold by default; (2) evbc = 2, wvbc = 3; or (3) periodic: ewpbc =\n\
    0. This is written in Matlab as:\n\
        c.dipdeg = -12; c.depth_min = 0; c.nvbc = 1; c.ewpbc = 0;\n\
\n\
  bc_periodic_nlayers [scalar]: Number of periodic source image layers to use if\n\
    any periodic BC are specified. 0 is just the primary source patch, which is\n\
    not a bad approximation because the primary source patch is chosen to be the\n\
    nearest periodically-repeated source to a receiver. A number K larger than 0\n\
    uses K layers of images. I recommend 1, which is the default if left\n\
    unspecified.\n\
\n\
\n\
"

#define _dc3dm_compress_PrintHelp_text_ \
"compress\n\
\n\
  build_read_filename [string]: the .rmu, .ra, .bmb, .hd, and .build files,\n\
    without suffix, from dc3dmBuild.\n\
\n\
  hm_write_filename [string]: The file to which to write the H-Matrix data. The\n\
    suffix .hm is appended to it.\n\
      In addition, a file with suffix .bc is created. It stores 4*Nr binary\n\
    double-precision numbers containing the velocity boundary condition\n\
    data. The BC is that the medium on all sides specified by the '[e|n|w|s]vbc'\n\
    data slide at the same speed in the direction of 'component'. If a side does\n\
    not have a velocity BC, then the entries are 0. Ordering is as follows: Nr\n\
    double-precisions for the E BC, then N, then W, then S.\n\
      Finally, a file with suffix .compress is created. It is a small key-value\n\
    file containing metadata.\n\
\n\
  hm_use_filename [string, optional]: An old H-matrix file for the same problem\n\
    that can be used to speed up constructing this one. This is generally useful\n\
    only if the old H-matrix was constructed at a higher tolerance than this new\n\
    one. Warning: dc3dmCompress does not check that all Green's functions\n\
    parameters are the same, so be careful when using this option.\n\
\n\
  allow_overwrite [scalar, optional]: Says whether the files associated with\n\
    hm_write_filename can be overwritten. If it is not provided, then overwrite\n\
    of these files is not allowed. If only the .bc file exists, then overwrite\n\
    is allowed.\n\
\n\
  src_disl [3-element vector]: Source dislocation in fault coordinates. This\n\
    vector corresponds to [DISL1, DISL2, DISL3] = [strike, dip, tensile]\n\
    dislocations in the notation of Okada's DC3D.\n\
\n\
  rcv_traction [3-element vector]: Compute tractions along this vector. The\n\
    vector is in fault coordinates. The components are [strike, dip, normal].\n\
    For example, to get along-dip tractions, set rcv_tractions = [0, 1, 0].\n\
\n\
  tol [scalar, optional]: Specifies the approximation error:\n\
        ||B_approx - B||_F <= tol ||B||_F,\n\
  where ||.||_F is the Frobenius norm, B is the matrix of Green's functions, and\n\
  B_approx is the H-matrix approximation. Default is 1e-5.\n\
\n\
  B_frobenius [scalar, optional]: Specifies ||B||_F, most likely an estimate. If\n\
    not provided, it is estimated. It is essential to provide B_frobenius when\n\
    computing off-diagonal blocks of a multi-component BEM matrix.\n\
\n\
  mu, nu [scalar]: Lame parameters.\n\
\n\
  nthreads [scalar, optional]: Number of (OpenMP) threads to use [default 1].\n\
\n\
  estimate_B_frobenius_only [0 or 1, optional]: If 1, estimate ||B||_F and\n\
    exit. The estimate is in the .compress file.\n\
\n\
\n\
"

#define _dc3dm_compressxkx_PrintHelp_text_ \
"compresskxk\n\
\n\
  Not yet.\n\
\n\
\n\
"

#endif
