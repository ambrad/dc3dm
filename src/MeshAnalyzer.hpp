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
 
#ifndef INCLUDE_UTIL_MESHANALYZER
#define INCLUDE_UTIL_MESHANALYZER

#include "util/include/OpenMP.hpp"
#include "RectMeshUnstruct_pri.hpp"
#include "Triangulation.hpp"

// Use Euclidean (2-norm) distance, rather than Manhattan (1-norm) distance,
// when computing neighborhoods and related quantities.
#define EUCLIDEAN_DISTANCE

//todo-remove when no longer interesting.
// Try out some ideas to bring src=srike, rcv=dip (and vice versa) problems up
// to OOA 2.
//#define SHEAR_SD_FACTOR

namespace util {
namespace rmesh {
using namespace std;

class MeshAnalyzer {
public:
  MeshAnalyzer(const RectMeshUnstruct* rmu, const RmuAnalyzer* ra);
  void ResetUnitFactor(size_t uf);

  const Rect& GetDomain () const { return _rmu->GetDomain(); }
  const Boundaries& GetBoundaries () const { return _ra->GetBoundaries(); }
  const vector<Rect>& GetRects () const { return _rmu->GetRects(); }
  int GetRectId (double x, double y) const { return _rmu->GetRectId(x, y); }

  const vector<RectId>& GetNbrs (RectId id) const { return _ra->GetNbrs(id); }
  const vector<Dir::Enum>& GetNbrDirs (RectId id) const
  { return _ra->GetNbrsBdys(id); }
  Dir::Enum GetBoundaryDir (RectId id) const
  { return _ra->GetBoundaryDir(id); }
  void GetCornerBC (Dir::Enum corner, Dir::Enum& bdy, Boundaries::BC& bc) const
  { return _ra->GetCornerBC(corner, bdy, bc); }
  void Distance (RectId id1, RectId id2, double* dx, double* dy,
                 Dir::Enum* dir = NULL) const
  { _ra->Distance(id1, id2, dx, dy, dir); }

  const Matrix<double>& GetX () const { return _ra->GetX(); }
  const Matrix<RectId>& GetTri () const { return _ra->GetTri(); }
  bool GetTriTag (double x, double y, TriTag* tag) const
  { return _ra->GetTriTag(x, y, tag); }
  void GetTriVertices (const TriTag& tag, double* v1, double* v2,
                       double* v3) const
  { _ra->GetPeriodicTri2(tag, v1, v2, v3); }

  size_t nr () const { return _nr; }
  size_t nx () const { return _nx; }
  // Use units of the smallest cell size.
  const vector<size_t>& Nc () const { return _nc; }
  // Unit cell size.
  double ux () const { return _ux; }
  double uy () const { return _uy; }
  // Convert to unit-cell units.
  int UnitDx(double dx) const;
  int UnitDy(double dy) const;
  double Distance(RectId id1, RectId id2, Dir::Enum* dir = NULL) const;
  double Distance(const Rect& r1, const Rect& r2) const;
  Rect& TranslateRect(Rect& r, Dir::Enum dir) const;
  Rect& TranslateRect(Rect& r, int Dx, int Dy) const;

  //testing
  const Matrix<Dir::Enum>& GetTriBdy () const { return _ra->GetTriBdy(); }

private:
  const RectMeshUnstruct* _rmu;
  const RmuAnalyzer* _ra;
  size_t _nr, _nx;
  // Unit cell size.
  double _ux, _uy;
  // Each cell can be expressed as a K x K square, where K is a multiple of
  // the unit cell size. _nc[id] is K for rect id.
  vector<size_t> _nc;

private:
  void Init(size_t uf = 1);
};

class InterpolatorMatrix {
public:
  InterpolatorMatrix(const MeshAnalyzer* ma) : _ma(*ma) {}
  virtual ~InterpolatorMatrix() {}

  // Number of layers of triangles used to interpolate a point in the center
  // triangle. The center triangle is layer 1. The next set of triangles
  // outward make layer 2, and so on.
  virtual size_t GetNbrTriLayersInSupport() const = 0;

  virtual void Call
  (// The function to interpolate, z, is 1 at node id_z1 and 0 everywhere else.
   size_t id_z1,
   // Interpolate to these points.
   const vector<double>& xi, const vector<double>& yi,
   // Point i is in tri(:, in[i].id).
   const vector<TriTag>& in,
   // Set zi[i] to the interpolated value at (xi[i], yi[i]).
   vector<double>& zi) const = 0;

  // For testing.
  virtual void GetSupportingVertices(size_t ti, vector<RectId>& vis) const = 0;

protected:
  const MeshAnalyzer& _ma;
};

struct LinearInterpolatorMatrix : public InterpolatorMatrix {
  LinearInterpolatorMatrix(const MeshAnalyzer* ma) : InterpolatorMatrix(ma) {}

  virtual size_t GetNbrTriLayersInSupport () const { return 1; };
  virtual void Call(size_t id_z1, const vector<double>& xi,
                    const vector<double>& yi, const vector<TriTag>& in,
                    vector<double>& zi) const;
  virtual void GetSupportingVertices(size_t ti, vector<RectId>& vis) const;
};

class CubicInterpolatorMatrix : public InterpolatorMatrix {
public:
  CubicInterpolatorMatrix(const MeshAnalyzer* ma);
  ~CubicInterpolatorMatrix();

  virtual size_t GetNbrTriLayersInSupport () const { return 2; }
  virtual void Call(size_t id_z1, const vector<double>& xi,
                    const vector<double>& yi, const vector<TriTag>& in,
                    vector<double>& zi) const;
  virtual void GetSupportingVertices(size_t ti, vector<RectId>& vis) const;

private:
  TriAnalyzer _ta;
  mutable vector<CiSupport*> _ss;
  mutable vector<omp_lock_t> _ss_lock;
  mutable vector<bool> _ss_set;

  void InitSupportIfNeeded(size_t ti) const;
};

}}

#include "MeshAnalyzer_inl.hpp"

#endif
