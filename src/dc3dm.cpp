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
 
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <sstream>
#include "util/include/Util.hpp"
#include "util/include/KeyValueFile.hpp"
#include "util/include/OpenMP.hpp"
#include "util/include/CodeAnalysis.hpp"
#include "hmmvp/include/Hd.hpp"
#include "hmmvp/include/Compress.hpp"
#include "Elastostatics.hpp"
#include "RectMeshUnstruct_pri.hpp"
#include "BrickMeshBemBuilder_pri.hpp"
#include "dc3dm.hpp"
using namespace util;
using namespace std;

//#define REPORT_MEMORY
#ifdef REPORT_MEMORY
#include <sys/resource.h>
#endif

#define fail_if_not(cond, msg) do {             \
    if (!(cond)) { err = msg; return false; }   \
  } while (0)

#define kvf_have_or_fail(fn, fld, out) do {                             \
    if (!kvf->fn(fld, out)) { err = "Missing " fld "."; return false; } \
  } while (0)

#define run_and_exit_if_Exception(code, announce)                       \
  try {                                                                 \
    cout << "-> " << announce << endl;                                  \
    code                                                                \
  } catch (const Exception& e) {                                        \
    cerr << "When " << announce << ": " << e.GetMsg() << endl;          \
    break;                                                              \
  }

static bool TestFileWrite (const string& fn) {
  FILE* fid = fopen(fn.c_str(), "wb");
  bool ret = fid;
  if (ret) fclose(fid);
  return ret;
}

namespace dc3dm {
bool Inputs::GetInputs (const string& kvf_fn) {
  KeyValueFile* kvf = NewKeyValueFile();
  if (!kvf->Read(kvf_fn) && !kvf->Read(kvf_fn + ".kvf")) {
    cerr << "Can't read " << kvf_fn << endl;
    DeleteKeyValueFile(kvf);
    return false;
  }
  bool ret = GetInputs(kvf);
  DeleteKeyValueFile(kvf);
  return ret;
}

bool Inputs::GetInputs (const KeyValueFile* kvf) {
  string err;
  if (!ProcessKvf(kvf, err)) {
    cerr << "Error: " << err << endl;
    return false;
  }
  return true;
}

namespace mesh {
void PrintHelp () {
  fprintf(stdout, _dc3dm_mesh_PrintHelp_text_);
  fprintf(stdout, _dc3dm_Header_text_);
}

struct Inputs : public dc3dm::Inputs {
  double min_len, max_len;
  Matrix<double> x, y;
  Matrix<double> f;
  string mesh_write_filename;
  size_t refine;
  bool do_uniform;

  virtual bool ProcessKvf (const KeyValueFile* kvf, string& err);
};

bool Ascending (const Matrix<double>& A) {
  for (size_t i = 2; i <= A.Size(); i++)
    if (A(i) <= A(i-1)) return false;
  return true;
}

bool AllPositive (const Matrix<double>& f, double f_min) {
  if (f_min > 0) return true;
  for (size_t i = 1; i <= f.Size(); i++)
    if (f(i) <= 0) return false;
  return true;
}

bool AnyNan (const Matrix<double>& f) {
  for (size_t i = 1; i <= f.Size(); i++)
    if (isnan(f(i))) return true;
  return false;
}

bool Inputs::ProcessKvf (const KeyValueFile* kvf, string& err) {
  Inputs& in = *this;
  const Matrix<double>* m;
  const string* s;
  double d;

  kvf_have_or_fail(GetMatd, "x", m);
  in.x = *m;
  fail_if_not(Ascending(in.x), "x must be monotonically increasing.");
  kvf_have_or_fail(GetMatd, "y", m);
  in.y = *m;
  fail_if_not(Ascending(in.y), "y must be monotonically increasing.");
  fail_if_not(in.x.Size() >= 2 && in.y.Size() >= 2,
              "numel(x), numel(y) >= 2");
  in.do_uniform = false;
  if (kvf->GetDouble("do_uniform", d))
    in.do_uniform = (bool) d;
  if (!in.do_uniform) {
    kvf_have_or_fail(GetMatd, "f", m);
    in.f = *m;
    fail_if_not(in.x.Size() == in.f.Size(2) && in.y.Size() == in.f.Size(1),
                "f must be numel(x) by numel(y).");
  }

  kvf_have_or_fail(GetDouble, "max_len", in.max_len);
  fail_if_not(in.max_len > 0,
              "max_len must be > 0 everywhere.");
  if (in.do_uniform)
    in.min_len = 0;
  else {
    kvf_have_or_fail(GetDouble, "min_len", in.min_len);
    fail_if_not(AllPositive(in.f, in.min_len),
                "max(f, min_len) must be > 0 everywhere.");
    fail_if_not(!AnyNan(in.f), "f contains NaN.");
  }
  fail_if_not(in.min_len <= in.max_len, "min_len must be <= max_len.");
  
  kvf_have_or_fail(GetString, "mesh_write_filename", s);
  in.mesh_write_filename = *s;
  in.refine = 0;
  if (kvf->GetDouble("refine", d) && d >= 0) in.refine = (size_t) d;
  return true;
}

void PrintProblemSetup (const rmesh::RectMeshUnstruct* rmu) {
  const vector<rmesh::Rect>& rs = rmu->GetRects();
  if (rs.empty()) {
    cout << "The mesh is empty. This should not happen." << endl;
    return;
  }
  double xs, xe, ys, ye, min_dx, max_dx, min_dy, max_dy;
  for (size_t i = 0; i < rs.size(); i++) {
    const rmesh::Rect& r = rs[i];
    if (i == 0) {
      xs = r.x; xe = r.x + r.dx;
      ys = r.y; ye = r.y + r.dy;
      min_dx = max_dx = r.dx;
      min_dy = max_dy = r.dy;
    } else {
      xs = min(xs, r.x); xe = max(xe, r.x + r.dx);
      ys = min(ys, r.y); ye = max(ye, r.y + r.dy);
      min_dx = min(min_dx, r.dx); max_dx = max(max_dx, r.dx);
      min_dy = min(min_dy, r.dy); max_dy = max(max_dy, r.dy);
    }
  }
  cout << "The domain is x: " << xs << " to " << xe << "; eta: " << ys << " to "
       << ye << ";" << endl
       << "  with " << rs.size() << " elements;" << endl
       << "  with min dx   " << min_dx << ", max dx   " << max_dx/min_dx
       << " * min dx" << endl
       << "   and min deta " << min_dy << ", max deta " << max_dy/min_dy
       << " * min deta;" << endl
       << "  with element aspect ratio " << max_dx/max_dy << "." << endl;
}

int Run (Inputs& in) {
  // Make sure we can write the output before doing any calculations.
  if (!TestFileWrite(in.mesh_write_filename + ".rmu")) {
    cerr << "Can't write " << in.mesh_write_filename + ".rmu" << endl;
    return -1;
  }

  const double xlen = in.x(in.x.Size()) - in.x(1);
  const double ylen = in.y(in.y.Size()) - in.y(1);
  const rmesh::Rect domain(in.x(1), in.y(1), xlen, ylen);
  // Represent the user's resolution function by a tensor-mesh discretization of
  // the domain with function values on the mesh points.
  rmesh::TensorMeshLinInterpRF* rf = NULL;

  if (in.do_uniform) {
    double
      x1 = in.x(1), xf = in.x(in.x.Size()),
      y1 = in.y(1), yf = in.y(in.y.Size());
    in.x.Resize(2); in.x(1) = x1; in.x(2) = xf;
    in.y.Resize(2); in.y(1) = y1; in.y(2) = yf;
    in.f.Resize(2, 2);
    for (size_t i = 1; i <= 4; i++) in.f(i) = in.max_len;
  }
  rf = new rmesh::TensorMeshLinInterpRF(in.x, in.y, in.f);

  // Options for meshing. Right now, just put lower and upper bounds on the
  // element size.
  const rmesh::RectOpts ro(in.min_len, in.max_len, in.do_uniform);

  rmesh::RectMeshUnstruct* rmu = NULL;
  do {
    // Discretize the domain.
    run_and_exit_if_Exception
      (rmu = NewRectMeshUnstruct(domain, ro, rf);,
       "constructing the RectMeshUnstruct");

    // For convergence testing, I need to refine the mesh. In fact, I do it
    // after smoothing in dc3dmBuild, but I'll leave this here for now.
    for (size_t i = 0; i < in.refine; i++) rmu->Refine();

    // Save the mesh to a file.
    run_and_exit_if_Exception
      (rmu->Serialize(in.mesh_write_filename + ".rmu");,
       "writing the RectMeshUnstruct (" + in.mesh_write_filename + ".rmu)");

    PrintProblemSetup(rmu);
  } while (0);

  if (rmu) rmesh::DeleteRectMeshUnstruct(rmu);
  if (rf) delete rf;

  return 0;
}

int Run (const string& kvf_fn) {
  // Read the key-value file.
  Inputs in;
  if (!in.GetInputs(kvf_fn)) return -1;
  return Run(in);
}
}

namespace build {
using namespace rmesh;

void PrintHelp () {
  fprintf(stdout, _dc3dm_build_PrintHelp_text_);
  fprintf(stdout, _dc3dm_Header_text_);
}

struct Inputs;

class BcProcessor {
public:
  void ProcessKvf(const KeyValueFile* kvf) throw (Exception);
  void ConstructBcs(const Inputs& in, const RectMeshUnstruct* rmu)
    throw (Exception);
  const Boundaries& GetBoundaries () const { return _b; }

private:
  Boundaries _b;

  struct Entry {
    double d;
    Dir::Enum dir;
    bool operator<(const Entry& e2) const { return d > e2.d; }
  };
  vector<Entry> _dominance;
};

struct Inputs : public dc3dm::Inputs {
  string mesh_read_filename, build_write_filename;
  double depth_min, dipdeg, neighborhood;
  size_t bc_periodic_nlayers;
  BcProcessor bcp;
  bool do_fullspace;
  size_t interp_type, refine;

  virtual bool ProcessKvf (const KeyValueFile* kvf, string& err);
};

Dir::Enum CharToDir (char d) {
  switch (d) {
  case 'e': return Dir::E;
  case 'w': return Dir::W;
  case 'n': return Dir::N;
  case 's': return Dir::S;
  default: return Dir::inside;
  }
}

ostream& operator<< (ostream& os, const Rect& r) {
  os << "[x: " << r.x << " to " << r.x + r.dx << ", eta: " << r.y << " to "
     << r.y + r.dy  << "]";
  return os;
}

string DirToString (Dir::Enum dir) {
  switch (dir) {
  case Dir::E: return " East";
  case Dir::N: return "North";
  case Dir::W: return " West";
  case Dir::S: return "South";
  default: assert(false); return "oops";
  }
}

string BCToString (Boundaries::BC bc) {
  switch (bc) {
  case Boundaries::bc_velocity:  return "    Velocity";
  case Boundaries::bc_0velocity: return "  0-Velocity";
  case Boundaries::bc_periodic:  return "    Periodic";
  case Boundaries::bc_free:      return "Free Surface";
  }
}

string ComponentToString (size_t c) {
  switch (c) {
  case 0: return "along-strike";
  case 1: return "along-dip";
  case 2: return "fault-normal";
  default: assert(false); return "oops";
  }
}

void BcProcessor::ProcessKvf (const KeyValueFile* kvf) throw (Exception) {
  const char* bc_strs[] = {"e0vbc", "n0vbc", "w0vbc", "s0vbc",
                           "evbc", "nvbc", "wvbc", "svbc",
                           "ewpbc", "nspbc", "wepbc", "snpbc"};
  const size_t n_bc_strs = sizeof(bc_strs) / sizeof(char*);
  _dominance.resize(4);
  for (size_t i = 0; i < 4; i++) _dominance[i].dir = (Dir::Enum) i;

  // Process the BC specification.
  for (size_t i = 0; i < n_bc_strs; i++) {
    double d;
    if (kvf->GetDouble(bc_strs[i], d)) {
      if (bc_strs[i][1] == '0' || bc_strs[i][1] == 'v') {
        Dir::Enum dir = CharToDir(bc_strs[i][0]);
        if (bc_strs[i][1] == '0') _b.SetZeroVelocityBC(dir);
        else {
          // Rect() is a stand-in for what comes later.
          _b.SetVelocityBC(dir, Rect());
        }
        if (d < 0) d = 0;
        _dominance[(int) dir].d = d;
        _dominance[(int) dir].dir = dir;
      } else {
        Dir::Enum dir = bc_strs[i][0] == 'e' || bc_strs[i][1] == 'e'
          ? Dir::E : Dir::N;
        _b.SetPeriodicBC(dir);
        _dominance[(int) dir].d = -1;
        _dominance[(int) dir].dir = dir;
        dir = Dir::Opposite(dir);
        _dominance[(int) dir].d = -1;
        _dominance[(int) dir].dir = dir;
      }
    }
  }
  std::sort(_dominance.begin(), _dominance.end());
}

inline double sind(double a) { return sin(a*M_PI/180.0); }

void BcProcessor::ConstructBcs (const Inputs& in, const RectMeshUnstruct* rmu)
  throw (Exception) {
  // Set up the v-BC slabs.
  const Rect& dmn = rmu->GetDomain();
  const double L = std::max(100.0, 10.0*(1 + in.bc_periodic_nlayers)) *
    std::max(dmn.dx, dmn.dy);
  // Prevent a v-BC slab from going the the surface.
  Dir::Enum surf_dir = in.dipdeg >= 0 ? Dir::N : Dir::S;
  double L_max = in.depth_min / fabs(sind(in.dipdeg));
  if (L_max > L) L_max = L;
  vector<bool> filled(9, false);
  for (size_t i = 0; i < 4; i++) {
    Dir::Enum dir = _dominance[i].dir;
    Dir::Enum ccw = Dir::OneCcw(dir), cw = Dir::OneCw(dir);
    if (_b.GetBC(dir) == Boundaries::bc_velocity) {
      const bool cw_filled = filled[(int) cw], ccw_filled = filled[(int) ccw];
      Rect r;
      switch (dir) {
      case Dir::E: {
        r.x = dmn.x + dmn.dx;
        r.dx = L;
        r.y = dmn.y - (cw_filled ? 0 : (surf_dir == Dir::S ? L_max : L));
        double y1 = dmn.y + dmn.dy +
          (ccw_filled ? 0 : (surf_dir == Dir::N ? L_max : L));
        r.dy = y1 - r.y;
      } break;
      case Dir::N: {
        r.y = dmn.y + dmn.dy;
        r.dy = (surf_dir == Dir::N ? L_max : L);
        r.x = dmn.x - (ccw_filled ? 0 : L);
        double x1 = dmn.x + dmn.dx + (cw_filled ? 0 : L);
        r.dx = x1 - r.x;
      } break;
      case Dir::W: {
        r.x = dmn.x - L;
        r.dx = L;
        r.y = dmn.y - (ccw_filled ? 0 : (surf_dir == Dir::S ? L_max : L));
        double y1 = dmn.y + dmn.dy +
          (cw_filled ? 0 : (surf_dir == Dir::N ? L_max : L));
        r.dy = y1 - r.y;
      } break;
      case Dir::S: {
        r.y = dmn.y - (surf_dir == Dir::S ? L_max : L);
        r.dy = surf_dir == Dir::S ? L_max : L;
        r.x = dmn.x - (cw_filled ? 0 : L);
        double x1 = dmn.x + dmn.dx + (ccw_filled ? 0 : L);
        r.dx = x1 - r.x;
      } break;
      default: assert(false);
      }
      _b.SetVelocityBC(dir, r);
    }
    if (Boundaries::IsVBC(_b.GetBC(dir)))
      filled[(int) dir] = filled[(int) cw] = filled[(int) ccw] = true;
  }

  // Check for a free boundary condition.
  if (in.depth_min == 0) _b.SetFreeBC(surf_dir);
}

hmmvp::Hd* ConstructHd (const MeshAnalyzer* ma) {
  Matrix<double> pbs(2, 3);
  pbs.Zero();
  const Rect& d = ma->GetDomain();
  if (ma->GetBoundaries().GetBC(Dir::E) ==
      Boundaries::bc_periodic) {
    pbs(1, 1) = d.x;
    pbs(2, 1) = d.x + d.dx;
  }
  if (ma->GetBoundaries().GetBC(Dir::N) ==
      Boundaries::bc_periodic) {
    pbs(1, 2) = d.y;
    pbs(2, 2) = d.y + d.dy;
  }

  // Rect centers in (x, eta) space. Could do it in (x, y, z) space, but it's
  // easier to specifiy the periodic boundaries in (x, eta) space.
  const vector<Rect>& rs = ma->GetRects();
  const Matrix<double>& xs = ma->GetX();
  Matrix<double> ctrs(3, rs.size());
  for (size_t i = 1; i <= rs.size(); i++) {
    ctrs(1, i) = xs(1, i);
    ctrs(2, i) = xs(2, i);
    ctrs(3, i) = 0;
  }

#ifdef TESTING_AND_ANALYSIS
  hmmvp::Hd* hd = hmmvp::NewHdAxisAligned(ctrs, &pbs, 2);
  printf("!!! Hd axis aligned !!!\n");
#else
  hmmvp::Hd* hd = hmmvp::NewHd(ctrs, &pbs);
#endif
  return hd;
}

// Use this as a place holder to serialize the BrickMeshBemBuilder.
struct UnusedGreensFn : public GreensFn {
  virtual double Call (const Rect&, double, double) const { return 0; }
};

bool Inputs::ProcessKvf (const KeyValueFile* kvf, string& err) {
  Inputs& in = *this;
  const string* s;
  double d;

  kvf_have_or_fail(GetString, "mesh_read_filename", s);
  in.mesh_read_filename = *s;
  kvf_have_or_fail(GetString, "build_write_filename", s);
  in.build_write_filename = *s;

  in.do_fullspace = false;
  if (kvf->GetDouble("do_fullspace", d)) in.do_fullspace = (bool) d;
  
  if (in.do_fullspace) {
    in.depth_min = 1; // Just need a nonzero value.
    in.dipdeg = 0;    // A horizontal fault is convenient.
  } else {
    kvf_have_or_fail(GetDouble, "depth_min", in.depth_min);
    fail_if_not(in.depth_min >= 0,
                "depth_min must be >= 0 (depth is a nonnegative number).");
    kvf_have_or_fail(GetDouble, "dipdeg", in.dipdeg);
  }

  in.neighborhood = 4;
  if (kvf->GetDouble("neighborhood", d)) {
    fail_if_not(in.neighborhood >= 0, "neighborhood must be >= 0.");
    in.neighborhood = d;
  }

  in.bc_periodic_nlayers = 1;
  if (kvf->GetDouble("bc_periodic_nlayers", d)) {
    fail_if_not(in.bc_periodic_nlayers >= 0 && in.bc_periodic_nlayers <= 50,
                "bc_periodic_nlayers must be >= 0 and <= 50.");
    in.bc_periodic_nlayers = (size_t) d;
  }
  try {
    in.bcp.ProcessKvf(kvf);
    fail_if_not(in.do_fullspace ||
                in.bcp.GetBoundaries().GetBC(Dir::N) != Boundaries::bc_periodic,
                "A N-S boundary condition can be specified only if the fault "
                "is in a full space.");
  } catch (const Exception& e) {
    err = e.GetMsg();
    return false;
  }

  in.interp_type = 3;
  if (kvf->GetDouble("interp_type", d)) {
    in.interp_type = (size_t) d;
    fail_if_not(in.interp_type == 1 || in.interp_type == 3,
                "interp_type must be 1 (linear) or 3 (cubic).");
  }

  in.refine = 0;
  if (kvf->GetDouble("refine", d) && d >= 0) in.refine = (size_t) d;

  return true;
}

void WriteBuildFile (const Inputs& in) throw (FileException) {
  KeyValueFile* kvf = NewKeyValueFile();
  kvf->AddDouble("depth_min", in.depth_min);
  kvf->AddDouble("dipdeg", in.dipdeg);
  kvf->AddDouble("do_fullspace", (double) in.do_fullspace);
  kvf->AddDouble("interp_type", (double) in.interp_type);
  if (!kvf->Write(in.build_write_filename + ".build"))
    throw FileException("Can't write " + in.build_write_filename + ".build");
  DeleteKeyValueFile(kvf);
}

void WriteElemFile (const vector<Rect>& rs, const string& bfn)
  throw (FileException)
{
  FILE* fid = fopen((bfn + ".elem").c_str(), "wa");
  if (!fid) throw FileException("Can't write " + bfn + ".elem");
  for (size_t i = 0; i < rs.size(); i++) {
    const Rect& r = rs[i];
    fprintf(fid, "%ld,%1.18e,%1.18e,%1.18e,%1.18e\n", i + 1, r.x + 0.5*r.dx,
            r.y + 0.5*r.dy, r.dx, r.dy);
  }
  fclose(fid);
}

void PrintProblemSetup (const Inputs& in, const MeshAnalyzer* ma) {
  cout << "Domain: " << ma->GetDomain() << "." << endl
       << "Boundary conditions:" << endl;
  bool have_periodic = false;
  for (size_t i = 0; i < 4; i++) {
    Dir::Enum dir = (Dir::Enum) i;
    Boundaries::BC bc = ma->GetBoundaries().GetBC(dir);
    if (bc == Boundaries::bc_periodic) have_periodic = true;
    cout << "  " << DirToString(dir) << ": " << BCToString(bc);
    if (bc == Boundaries::bc_velocity)
      cout << " with slab " << ma->GetBoundaries().GetRect(dir);
    cout << "." << endl;
  }
  if (have_periodic)
    cout << "Using " << in.bc_periodic_nlayers << " layer"
         << (in.bc_periodic_nlayers != 1 ? "s" : "")
         << " of periodic source images." << endl;
  Dir::Enum surf_dir = in.dipdeg >= 0.0 ? Dir::N : Dir::S;
  if (in.do_fullspace)
    cout << "Fullspace." << endl;
  else {
    cout << "Halfspace. The updip side is " << DirToString(surf_dir)
         << "; it is at depth " << in.depth_min << "." << endl;
    cout << "The fault dips at " << fabs(in.dipdeg) << " degrees." << endl;
  }
  cout << "The neighborhood factor is " << in.neighborhood << "." << endl;
}

int Run (Inputs& in) {
  // Make sure we can write the output before doing any calculations.
  if (!TestFileWrite(in.build_write_filename + ".ra")) {
    cerr << "Can't write " << in.build_write_filename + ".*" << endl;
    return -1;
  }

  rmesh::RectMeshUnstruct* rmu = NULL;
  rmesh::RectMeshUnstruct* srmu = NULL;
  rmesh::RmuAnalyzer* ra = NULL;
  rmesh::MeshAnalyzer* ma = NULL;
  hmmvp::Hd* hd = NULL;
  rmesh::InterpolatorMatrix* im = NULL;
  UnusedGreensFn unused;
  rmesh::BrickMeshBemBuilder* bmb = NULL;

  do {
    // Read in the mesh that dc3dmMesh made.
    run_and_exit_if_Exception
      (rmu = rmesh::NewRectMeshUnstruct(in.mesh_read_filename + ".rmu");,
       "reading the RectMeshUnstruct (" + in.mesh_read_filename + ".rmu)");

    // Smooth the mesh.
    run_and_exit_if_Exception
      (srmu = rmesh::SmoothRectMeshUnstruct(rmu, in.bcp.GetBoundaries());,
       "smoothing the RectMeshUnstruct");
    cout << "  (Went from " << rmu->GetRects().size() << " to "
         << srmu->GetRects().size() << " elems.)" << endl;
    rmesh::DeleteRectMeshUnstruct(rmu);
    rmu = NULL;

    // For convergence testing, I need to refine the mesh.
    if (in.refine > 0) {
      size_t no = srmu->GetRects().size();
      for (size_t i = 0; i < in.refine; i++) srmu->Refine();
      cout << "  (Refined from " << no << " to " << srmu->GetRects().size()
           << " elems." << endl;
    }

    run_and_exit_if_Exception
      (srmu->Serialize(in.build_write_filename + ".rmu");,
       "writing the smoothed RectMeshUnstruct (" + in.build_write_filename +
       ".rmu)");

    // Build the slabs that implement the velocity boundary conditions.
    run_and_exit_if_Exception
      (in.bcp.ConstructBcs(in, srmu);,
       "constructing the boundary conditions");

    // Make the mesh analyzer. The analyzer builds an adjacency graph for the
    // mesh elements and a triangulation based on element centers that is used
    // for interpolation.
    run_and_exit_if_Exception
      (ra = rmesh::NewRmuAnalyzer(srmu, in.bcp.GetBoundaries());,
       "constructing the RmuAnalyzer");

    // Write the RmuAnalyzer.
    run_and_exit_if_Exception
      (ra->Serialize(in.build_write_filename + ".ra");,
       "writing the RmuAnalyzer (" + in.build_write_filename + ".ra)");

    // The MeshAnalyzer is a wrapper container for the mesh objects.
    run_and_exit_if_Exception
      (ma = new rmesh::MeshAnalyzer(srmu, ra);,
       "constructing the MeshAnalyzer");

    // Build the H-matrix spatial decomposition Hd for this domain.
    run_and_exit_if_Exception
      (hd = ConstructHd(ma);,
       "constructing the Hd");

    // Write the Hd.
    run_and_exit_if_Exception
      (hmmvp::WriteHd(hd, in.build_write_filename + ".hd");,
       "writing the Hd (" + in.build_write_filename + ".hd)");

    // The InterpolationMatrix is used in the IGA method to compute elements of
    // G_IGA = A*G*I.
    run_and_exit_if_Exception
      (im = in.interp_type == 1 ?
       (rmesh::InterpolatorMatrix*) (new rmesh::LinearInterpolatorMatrix(ma)) :
       (rmesh::InterpolatorMatrix*) (new rmesh::CubicInterpolatorMatrix(ma));,
       "constructing the " + string(in.interp_type == 1 ? "Linear" : "Cubic") +
       " InterpolatorMatrix");

    // Make the BrickMeshBemBuilder. It uses the MeshAnalyzer to determine the
    // accuracy of quadrature to do for each source-receiver pair. dc3dmBuild
    // doesn't do any Green's function calculations, so we use an empty
    // implementation to get BrickMeshBemBuilder to do just the mesh analysis
    // part and quit.
    //   If the mesh is uniform, override the user's neighborhood because
    // there's no need to have it be > 0.
    if (srmu->GetRectOpts().uniform) in.neighborhood = 0;
    run_and_exit_if_Exception
      (UnusedGreensFn unused;
       bmb = rmesh::NewBrickMeshBemBuilder
       (&unused, ma, im, in.neighborhood, in.bc_periodic_nlayers);,
       "constructing the BrickMeshBemBuilder");

    // Write the BrickMeshBemBuilder setup data.
    run_and_exit_if_Exception
      (bmb->Serialize(in.build_write_filename + ".bmb");,
       "writing the BrickMeshBemBuilder (" + in.build_write_filename + ".bmb)");

    // Write a key-value file to carry some parameters over to dc3dmCompress.
    run_and_exit_if_Exception
      (WriteBuildFile(in);,
       "writing the build file (" + in.build_write_filename + ".build)");

    // Write the .elem file. Its purpose is to make the Matlab interface
    // optional. It provides the necessary mesh and element ordering
    // information.
    run_and_exit_if_Exception
      (WriteElemFile(srmu->GetRects(), in.build_write_filename);,
       "writing the element file (" + in.build_write_filename + ".elem)");

    PrintProblemSetup(in, ma);
  } while (0);

  if (bmb) DeleteBrickMeshBemBuilder(bmb);
  if (im) delete im;
  if (hd) hmmvp::DeleteHd(hd);
  if (ma) delete ma;
  if (ra) rmesh::DeleteRmuAnalyzer(ra);
  if (srmu) rmesh::DeleteRectMeshUnstruct(srmu);
  if (rmu) rmesh::DeleteRectMeshUnstruct(rmu);

  return 0;
}

int Run (const string& kvf_fn) {
  // Read the key-value file.
  Inputs in;
  if (!in.GetInputs(kvf_fn)) return -1;
  return Run(in);
}
}

namespace compress {
using namespace rmesh;

/* DEV A specific rectangle-to-point Green's function implements this interface
   for BrickMeshBemBuilder. If dc3dm is built to use OpenMP, then Call must be
   thread safe. */
struct ImplGreensFn : public rmesh::GreensFn {
  virtual ~ImplGreensFn() {}
  /* DEV Compute B(rs,cs). Indexing starts at 1. B is preallocated.
     Return true if all is well; false if there is an error in computing the
     Green's function and you want compression to stop. */
  virtual double Call(const rmesh::Rect& src, double rx, double ry) const = 0;
};

// Call Okada's routine dc3d for a rectangular dislocation in an elastic
// half-space.
class Dc3dGreensFn : public ImplGreensFn {
public:
  Dc3dGreensFn (
    const rmesh::Rect& domain, double depth_min, double dipdeg, double mu,
    double nu, bool do_fullspace, const double src_disl[3],
    const double rcv_traction[3])
    : _domain(domain), _lp(mu, nu), _dipdeg(dipdeg), _do_fullspace(do_fullspace)
  {
    _cd = es::cosd(_dipdeg);
    _sd = es::sind(_dipdeg);
    _depth0 = depth_min + fabs(_sd) * _domain.dy;
    _depth1 = depth_min;
    if (_sd < 0.0) std::swap(_depth0, _depth1);
    memcpy(&_src_disl[0], &src_disl[0], 3*sizeof(src_disl[0]));
    // Set up the vector along which to compute traction. We only need dip to be
    // right for the calculation to work.
    es::dc3::Elem es(0, _dipdeg, _cd, _sd, 1, 1, 1, 1, 0, 0);
    es.ToGlobal(rcv_traction, _rcv_along_global);
  }

  virtual double Call (const rmesh::Rect& src, double rx, double ry) const {
    Ca::GetTimer()->Tic(1);
    double eta0, alpha;
    eta0 = src.y - _domain.y;
    alpha = eta0 / _domain.dy;
    es::dc3::Elem es((1.0 - alpha) * _depth0 + alpha * _depth1, _dipdeg, _cd,
                     _sd, 0, src.dx, 0, src.dy, src.x, _cd * eta0);
    double s[6];
    { double obs[3], u[3], du[9];
      eta0 = ry - _domain.y;
      alpha = eta0 / _domain.dy;
      obs[0] = rx;
      obs[1] = eta0 * _cd;
      obs[2] = -((1.0 - alpha) * _depth0 + alpha * _depth1);
      es::dc3::Dc3d(_lp, es, _src_disl, obs, u, du, _do_fullspace);
      es::DuToS(_lp, du, s); }
    double B;
    es::ProjectStress(s, es.Normal(), _rcv_along_global, NULL, &B, NULL, NULL);
    Ca::GetTimer()->Toc(1);
    return B;
  }

private:
  rmesh::Rect _domain;
  es::LameParms _lp;
  double _depth0, _depth1, _dipdeg, _cd, _sd, _src_disl[3],
    _rcv_along_global[3];
  bool _do_fullspace;
  int _rcv_component;
};

void PrintHelp () {
  fprintf(stdout, _dc3dm_compress_PrintHelp_text_);
  fprintf(stdout, _dc3dm_Header_text_);
}

struct Inputs : public dc3dm::Inputs {
  string build_read_filename, hm_write_filename, hm_use_filename;
  bool allow_overwrite, do_fullspace, est_Bfro_only;
  double src_disl[3], rcv_traction[3];
  size_t nthreads, interp_type;
  double tol, mu, nu, depth_min, dipdeg;
  double Bfro;
  // Undoc'ed options for testing and analysis.
  bool use_bmb_lra, allow_0rank_blocks, use_compress_qr;
  hmmvp::Compressor::TolMethod tol_method;

  virtual bool ProcessKvf (const KeyValueFile* kvf, string& err);
};

bool CanRead(const string& fn) {
  FILE* fid = fopen(fn.c_str(), "rb");
  if (!fid) return false;
  fclose(fid);
  return true;
}

bool Inputs::ProcessKvf (const KeyValueFile* kvf, string& err) {
  Inputs& in = *this;
  const string* s;
  const Matrix<double>* m;
  double d, d1;

#define rfn(fn) kvf_have_or_fail(GetString, #fn, s); in.fn = *s;
  rfn(build_read_filename);
  rfn(hm_write_filename);
  in.hm_use_filename = "";
  if (kvf->GetString("hm_use_filename", s)) in.hm_use_filename = *s;
#undef rfn
  fail_if_not(in.hm_write_filename != in.hm_use_filename,
              "hm_write_filename can't be the same as hm_use_filename.");
  in.allow_overwrite = false;
  if (kvf->GetDouble("allow_overwrite", d)) in.allow_overwrite = (bool) d;
  fail_if_not(in.allow_overwrite || !CanRead(in.hm_write_filename + ".hm"),
              in.hm_write_filename + ".hm file exists; won't overwrite.");
  
  kvf_have_or_fail(GetMatd, "src_disl", m);
  fail_if_not(m->Size() == 3, "src_disl must be a 3-element vector.");
  memcpy(&in.src_disl[0], m->GetPtr(), 3*sizeof(double));
  kvf_have_or_fail(GetMatd, "rcv_traction", m);
  fail_if_not(m->Size() == 3, "rcv_traction must be a 3-element vector.");
  memcpy(&in.rcv_traction[0], m->GetPtr(), 3*sizeof(double));

  in.tol = 1e-5;
  if (kvf->GetDouble("tol", d)) {
    in.tol = d;
    fail_if_not(in.tol > 0 && in.tol <= 0.1, "tol must be > 0 and <= 0.1.");
  }

  in.Bfro = -1;
  if (kvf->GetDouble("B_frobenius", d)) {
    in.Bfro = d;
    fail_if_not(in.Bfro > 0, "B_frobenius must be > 0.");
  }
  
  kvf_have_or_fail(GetDouble, "mu", in.mu);
  fail_if_not(in.mu > 0, "mu must be > 0.");
  kvf_have_or_fail(GetDouble, "nu", in.nu);
  fail_if_not(in.nu >= 0, "nu must be >= 0.");
  fail_if_not(in.nu < 1, "nu must be < 1.");

  in.nthreads = 1;
  if (kvf->GetDouble("nthreads", d)) {
    if ((int) d <= 0) throw Exception("nthreads must be >= 1.");
    in.nthreads = (size_t) d;
    if (in.nthreads > 256) {
      cout << "Warning: nthreads > 256. That doesn't seem right. "
           << "Proceeding anyway." << endl;
    }
  }

  in.est_Bfro_only = false;
  if (kvf->GetDouble("estimate_B_frobenius_only", d))
    in.est_Bfro_only = (bool) d;

  KeyValueFile* kvf_build = NewKeyValueFile();
  fail_if_not(kvf_build->Read(in.build_read_filename + ".build"),
              "Can't read " + in.build_read_filename + ".build");
  fail_if_not(kvf_build->GetDouble("depth_min", in.depth_min) &&
              kvf_build->GetDouble("dipdeg", in.dipdeg) &&
              kvf_build->GetDouble("do_fullspace", d) &&
              kvf_build->GetDouble("interp_type", d1),
              "Failed when reading " + in.build_read_filename +
              ".build; should not happen.");
  in.do_fullspace = (bool) d;
  in.interp_type = (size_t) d1;
  DeleteKeyValueFile(kvf_build);

  // Undocumented keys for analysis and testing. These are the defaults and the
  // values I intend users to always use:
  in.use_bmb_lra = true;
  in.use_compress_qr = true;
  in.allow_0rank_blocks = false;
  in.tol_method = hmmvp::Compressor::tm_mrem_fro;
  // But for analysis and testing:
  if (kvf->GetDouble("use_bmb_lra", d)) in.use_bmb_lra = (bool) d;
  if (kvf->GetDouble("use_compress_qr", d)) in.use_compress_qr = (bool) d;
  if (kvf->GetDouble("allow_0rank_blocks", d)) in.allow_0rank_blocks = (bool) d;
  if (kvf->GetString("tol_method", s)) {
    if (*s == "brem-fro") in.tol_method = hmmvp::Compressor::tm_brem_fro;
    else if (*s == "mrem-abs") in.tol_method = hmmvp::Compressor::tm_mrem_abs;
    // Otherwise, stick with tm_mrem_fro.
  }

  return true;
}

static string Vec3ToStr (const double v[3]) {
  stringstream ss;
  ss << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
  return ss.str();
}

void PrintProblemSetup (const Inputs& in) {
  cout << "Reading from " << in.build_read_filename << " and "
       << in.build_read_filename << endl
       << "  and writing to " << in.hm_write_filename << "." << endl;
  if (in.hm_use_filename.empty())
    cout << "Not using an old H-matrix." << endl;
  else
    cout << "Using the old H-matrix " << in.hm_use_filename << "." << endl;
  omp_set_num_threads(in.nthreads);
  int nthreads = omp_get_max_threads();
  cout << "Using " << nthreads << " OpenMP threads (" << in.nthreads
       << " threads requested)." << endl
       << "H-matrix tolerance is ||G_true - G_approx||_F <= " << in.tol
       << " ||G_true||_F." << endl
       << "Source dislocation " << Vec3ToStr(in.src_disl).c_str() << "; "
       << "receiver traction component " << Vec3ToStr(in.rcv_traction).c_str()
       << "." << endl << "mu is " << in.mu << " and nu is " << in.nu << "."
       << endl;
  if (in.do_fullspace)
    cout << "  Fullspace." << endl;
  else
    cout << "  Halfspace with depth_min " << in.depth_min << " dip "
         << in.dipdeg << " degrees." << endl;
  if (!in.use_bmb_lra) cout << "Undoc'ed: not using bmb lra." << endl;
  if (!in.use_compress_qr) cout << "Undoc'ed: not using compress QR." << endl;
  if (in.allow_0rank_blocks)
    cout << "Undoc'ed: allowing rank-0 blocks." << endl;
  if (in.tol_method != hmmvp::Compressor::tm_mrem_fro)
    cout << "Undoc'ed: using tol_method " << in.tol_method << "." << endl;
}

// Compute the influence of velocity boundary condition large rectangles on the
// fault elements.
void ComputeBcs (const Inputs& in, BrickMeshBemBuilder* bmb, const size_t nr)
  throw (FileException) {
  Matrix<double> B(nr, 4);
  B.Zero();

  omp_set_num_threads(in.nthreads);
  for (size_t j = 0; j < 4; j++)
#pragma omp parallel for
    for (size_t i = 0; i < nr; i++)
      B(i+1, j+1) = bmb->Call(i, (Dir::Enum) j);
    
  FILE* fid = fopen((in.hm_write_filename + ".bc").c_str(), "wb");
  if (!fid)
    throw FileException("Can't write file " + in.hm_write_filename + ".bc");
  write(B.GetPtr(), 4*nr, fid);
  fclose(fid);
}

// Implement hmmvp::Compress's Green's function interface.
class HmmvpGreensFn : public hmmvp::GreensFn {
public:
  HmmvpGreensFn (BrickMeshBemBuilder* bmb) : _bmb(bmb), _alpha(0.1) {}

  virtual bool Call (const hmmvp::CompressBlockInfo& cbi,
                     const vector<UInt>& rs, const vector<UInt>& cs,
                     double* B) const {
    /* We want to tell BMB that it can make an error in computing this
       element. If each element is in error by E11, then the additional error in
       the block is
           norm(E, fro) = sqrt(m n) E11.
       We want
           sqrt(m n) E11 <= _alpha tol,
       and so
           E11 <= _alpha tol / sqrt(m n). */
    double tol = _alpha * cbi.tol / sqrt(cbi.nrows * cbi.ncols);
    for (size_t ic = 0, k = 0; ic < cs.size(); ic++)
      for (size_t ir = 0; ir < rs.size(); ir++, k++)
        B[k] = _bmb->Call(rs[ir] - 1, cs[ic] - 1, tol);
    return true;
  }

  virtual void StartBlock (const hmmvp::CompressBlockInfo& cbi) const {
    _bmb->StartBlock();
  }
  virtual void FinishBlock (const hmmvp::CompressBlockInfo& cbi) const {
    _bmb->FinishBlock();
  }

private:
  BrickMeshBemBuilder* _bmb;
  double _alpha;
};

void WriteCompressMetadata (const Inputs& in, const BrickMeshBemBuilder* bmb,
                            const hmmvp::Compressor* c, const string& fn) {
  KeyValueFile* kvf = NewKeyValueFile();
  kvf->AddDouble("B_frobenius_estimate", c->GetBfroEstimate());
  Matrix<double> m3(1, 3);
  for (size_t i = 0; i < 3; i++) m3(i+1) = in.src_disl[i];
  kvf->AddMatd("src_disl", m3);
  for (size_t i = 0; i < 3; i++) m3(i+1) = in.rcv_traction[i];
  kvf->AddMatd("rcv_traction", m3);
  kvf->AddDouble("fullspace", in.do_fullspace);
  kvf->AddDouble("interp_type", in.interp_type);
  kvf->AddDouble("tol", in.tol);
  kvf->AddDouble("mu", in.mu);
  kvf->AddDouble("nu", in.nu);
  kvf->AddDouble("Bfro_was_provided", in.Bfro > 0);
  bool ret = kvf->Write(fn);
  DeleteKeyValueFile(kvf);
  if (!ret) throw FileException("Can't write " + fn);
}

void PrintBmbStats (const BrickMeshBemBuilder* bmb, const MeshAnalyzer* ma) {
  size_t ns[2];
  ns[0] = bmb->GetCtrNcalls();
  ns[1] = bmb->GetCtrNgfs();

  size_t max_nc = 0;
  for (size_t i = 0; i < ma->Nc().size(); i++)
    max_nc = std::max(max_nc, ma->Nc()[i]);
  const Boundaries& b = ma->GetBoundaries();
  size_t pf = 0;
  if (b.GetBC(Dir::E) == Boundaries::bc_periodic) pf++;
  if (b.GetBC(Dir::N) == Boundaries::bc_periodic) pf++;
  switch (pf) {
  case 0: pf = 1; break;
  case 1: pf = 2 * bmb->GetBcPeriodicNlayers() + 1; break;
  case 2: pf = 2 * bmb->GetBcPeriodicNlayers() + 1; pf *= pf; break;
  default: assert (false); break;
  }
  printf("calls: std %lld IGA %lld (%1.2f)\n",
         (long long int) (pf * ns[0]), (long long int) ns[1],
         (double) ns[1] / (pf * ns[0]));
  printf("nrgfs: %ld\n", bmb->GetCtrNrgfs());
#ifdef ANALYZE_CODE
  double et_tot = Ca::GetTimer()->TotEt(0), et_call = Ca::GetTimer()->TotEt(1);
  printf("et tot %1.3f call %1.3f (%1.4f)\n", et_tot, et_call, et_call/et_tot);
#endif
}

int Run (const Inputs& in) {
  PrintProblemSetup(in);

  rmesh::RectMeshUnstruct* rmu = NULL;
  rmesh::RmuAnalyzer* ra = NULL;
  rmesh::MeshAnalyzer* ma = NULL;
  rmesh::InterpolatorMatrix* im = NULL;
  ImplGreensFn* igf = NULL;
  rmesh::BrickMeshBemBuilder* bmb = NULL;
  hmmvp::Hd* hd = NULL;
  HmmvpGreensFn* hgf = NULL;
  hmmvp::Compressor* c = NULL;

  do {
    // Read in the mesh that dc3dmMesh made.
    run_and_exit_if_Exception
      (rmu = rmesh::NewRectMeshUnstruct(in.build_read_filename + ".rmu");,
       "reading the RectMeshUnstruct (" + in.build_read_filename + ".rmu)");

    // Read in the .ra, .bmb, and .hd files from dc3dBuild and construct the
    // other dc3dBuild objects.
    run_and_exit_if_Exception
      (ra = NewRmuAnalyzer(rmu, in.build_read_filename + ".ra");,
       "reading the RmuAnalyzer (" + in.build_read_filename + ".ra)");
    run_and_exit_if_Exception
      (ma = new rmesh::MeshAnalyzer(rmu, ra);,
       "constructing the MeshAnalyzer");
    run_and_exit_if_Exception
      (im = in.interp_type == 1 ?
       (rmesh::InterpolatorMatrix*) (new rmesh::LinearInterpolatorMatrix(ma)) :
       (rmesh::InterpolatorMatrix*) (new rmesh::CubicInterpolatorMatrix(ma));,
       "constructing the " + string(in.interp_type == 1 ? "Linear" : "Cubic") +
       "InterpolatorMatrix");
    //todo If there are other Green's functions available, this part needs to be
    // made into a switch.
    run_and_exit_if_Exception
      (igf = new Dc3dGreensFn(
        ma->GetDomain(), in.depth_min, in.dipdeg, in.mu, in.nu, in.do_fullspace,
        in.src_disl, in.rcv_traction);,
       "constructing the BrickMeshBemBuilder Green's function");
    run_and_exit_if_Exception
      (bmb = rmesh::NewBrickMeshBemBuilder(igf, ma, im,
                                           in.build_read_filename + ".bmb");,
       "reading the BrickMeshBemBuilder (" + in.build_read_filename + ".bmb)");
    run_and_exit_if_Exception
      (hd = hmmvp::NewHd(in.build_read_filename + ".hd");,
       "reading the Hd (" + in.build_read_filename + ".hd)");

    hgf = new HmmvpGreensFn(bmb);

    run_and_exit_if_Exception
      (c = hmmvp::NewCompressor(hd, hgf);,
       "constructing the Compressor");
    c->SetOutputLevel(1);
    c->SetTolMethod(in.tol_method);
    c->SetTol(in.tol);
    c->Allow0RankBlocks(in.allow_0rank_blocks);
    bmb->SetUseLra(in.use_bmb_lra);
    bmb->SetOmpNthreads(in.nthreads);
    c->SetOmpNthreads(in.nthreads);
    c->UseCompressQr(in.use_compress_qr);
    if (!in.hm_use_filename.empty())
      run_and_exit_if_Exception
        (c->UseHmatFile(in.hm_use_filename + ".hm");,
         "using the old H-matrix " + in.hm_use_filename + ".hm");
    if (in.tol_method == hmmvp::Compressor::tm_mrem_fro) {
      if (in.Bfro > 0)
        c->SetBfroEstimate(in.Bfro);
      else {
        if (c->HaveOldHmat())
          c->SetBfroEstimate(c->GetOldHmatBfro());
        else {
          bmb->UsePoorQuality();
          run_and_exit_if_Exception
            (c->SetBfroEstimate(c->EstimateBfro());,
             "estimating ||B||_F");
          bmb->UseFullQuality();
        }
      }
    }
    run_and_exit_if_Exception
      (WriteCompressMetadata(in, bmb, c, in.hm_write_filename + ".compress");,
       "writing to " + in.hm_write_filename + ".compress");
    if (in.est_Bfro_only) {
      cout << "Exiting because estimate_B_frobenius_only." << endl;
      break;
    }
    run_and_exit_if_Exception
      (ComputeBcs(in, bmb, ma->nr());,
       "computing the boundary conditions and writing to " +
       in.hm_write_filename + ".bc");
    Ca::GetTimer()->Tic(0);
    run_and_exit_if_Exception
      (c->CompressToFile(in.hm_write_filename + ".hm");,
       "compressing to " + in.hm_write_filename + ".hm");
    Ca::GetTimer()->Toc(0);

    PrintBmbStats(bmb, ma);
  } while (0);

  if (c) hmmvp::DeleteCompressor(c);
  if (hgf) delete hgf;
  if (hd) hmmvp::DeleteHd(hd);
  if (bmb) DeleteBrickMeshBemBuilder(bmb);
  if (igf) delete igf;
  if (im) delete im;
  if (ma) delete ma;
  if (ra) rmesh::DeleteRmuAnalyzer(ra);
  if (rmu) rmesh::DeleteRectMeshUnstruct(rmu);

#ifdef REPORT_MEMORY
  rusage r;
  getrusage(RUSAGE_SELF, &r);
  errpr("mem %1.6f MB\n", r.ru_maxrss / (double) (2 << 10));
#endif

  return 0;
}

int Run (const string& kvf_fn) {
  // Read the key-value file.
  Inputs in;
  if (!in.GetInputs(kvf_fn)) return -1;
  return Run(in);
}
}

namespace compresskxk {
struct Inputs : public dc3dm::Inputs {
  virtual bool ProcessKvf (const util::KeyValueFile* kvf, std::string& err) {}
};

void PrintHelp () {
  fprintf(stdout, _dc3dm_compressxkx_PrintHelp_text_);
  fprintf(stdout, _dc3dm_Header_text_);
}

int Run (const Inputs& in) {
  cout << "Not yet; sorry." << endl;
  return 0;
}

int Run (const string& kvf_fn) {
  // Read the key-value file.
  Inputs in;
  if (!in.GetInputs(kvf_fn)) return -1;
  return Run(in);
}
}}
