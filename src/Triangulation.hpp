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
 
#ifndef INCLUDE_UTIL_TRIANGULATION
#define INCLUDE_UTIL_TRIANGULATION

#include "util/include/OpenMP.hpp"
#include "RectMeshUnstruct_pri.hpp"
#include "Tri2.hpp"

namespace util {
namespace rmesh {
using namespace std;

class MeshAnalyzer;

// Tri ids are base 1; rect (vertex) ids are base 0.
class TriAnalyzer {
public:
  TriAnalyzer (const MeshAnalyzer* ma) : _ma(*ma) { Init(); }

  // The edge is defined by (vertex vi1, vertex vi2). Return the triangle other
  // than ti that has this edge.
  size_t GetTriSharingEdge(size_t ti, RectId vi1, RectId vi2) const;
  // Edge i in 0:2 (with edge ordering) of triangle ti.
  size_t GetTriSharingEdge(size_t ti, size_t ei) const;

  const vector<size_t>& GetTrisSharingVertex(RectId vi) const
  { return _vi2ti[vi]; }
  
private:
  const MeshAnalyzer& _ma;
  vector< vector<size_t> > _vi2ti;

  void Init();
  void InitVertexToTrisMap();
};

/* A CiSupport contains a triangle, its (up to) three edge-sharing neighbors,
   its other vertex-sharing neighbors, and cubic interpolation data.
     Edge order is (vertex 0, 1), (1, 2), (2, 0).
     Vertex order is 0, 1, 2.
     Interpolation formulas come from Section 2 of [1] P. Alfeld, A trivariate
   Clough-Tocher scheme, 1984, which reviews in detail the classic bivariate
   Clough-Tocher splitting method. */
struct Connectivity;
class CiSupport {
public:
  CiSupport (const MeshAnalyzer* ma, const TriAnalyzer& ta, size_t ti);
  ~CiSupport();

  void GetSupportingVertices(size_t ti, vector<RectId>& vis) const
  { vis = _all_vis; }
  // The 1-function is 1 at one vertex and 0 at every other vertex. In this
  // case, it is 1 at vertex vi_is_one. Get the interpolated value at (xi[0],
  // xi[1]).
  double InterpWithOneFn(RectId vi_is_one, const TriTag& tag,
                         const double xi[2]);

private:
  const MeshAnalyzer& _ma;
  size_t _ti; // The primary triangle.

  // There are three vertex types to consider: those belonging to the primary
  // triangle (_ti); those belonging to the <= 3 triangles sharing an edge with
  // the primary triangle; and the rest. Each gives rise to a certain
  // coefficient pattern. Distinguish among these.
  vector<RectId> _pri_vis; // Primary triangle's vertices.
  vector<RectId> _es_vis;  // Those of the edge-sharing tris.
  // _all_vis = [_pri_vis _es_vis (the rest)]. This sets up a local vertex
  // numbering. We'll use lvi for local vertex index and vi for global. _all_vis
  // contains more than just the vertex nbrs as defined by the triangulations;
  // it contains all rect mesh vertex neighbors.
  vector<RectId> _all_vis;

  Matrix<size_t> _t; // Ordered [primary, edge-sharing, the rest], though we
                     // don't use this fact.
  Matrix<double> _x; // Ordered according to _all_vis.
  vector< vector<size_t> > _vi2ti; // Vertex-to-tri map indexing _t.
  // _lvi_ats[i] are the (primary tri) lvis affected by a 1-value at i.
  vector< vector<RectId> > _lvi_ats;

  // Barycentric data for the microtriangles in the primary triangle.
  double _x4[2]; // The primary triangle's centroid.
  Tri2 _tri[3];
  double _Ti[3][4];

  // The 19 coefficients associated with interpolating the 1-function. _cs[i]
  // are the coefficients when vertex _all_vis[i] is 1.
  vector<double*> _cs;
  vector<omp_lock_t> _cs_lock;
  vector<bool> _cs_set;

  // Ordering of coefficients is according to this enum:
  enum CoefIdx {
    c3000 = 0, c2100, c2010, c0300, c1200, c0210, c0030, c1020, c0120, c0021,
    c0201, c2001, c0111, c1011, c1101, c0012, c0102, c1002, c0003 };
  static const size_t _interp_coefs[19];
  static const size_t _lam_pow[19][4];
  static const size_t _nc = sizeof(_interp_coefs)/sizeof(_interp_coefs[0]);

  // Gradient estimate data. For fitting a quadratic, we need vertex adjacency.
  vector< vector<RectId> > _ag;

private:
  void SetConnectivity(const TriAnalyzer& ta, Connectivity& c);
  void SetTris(const Connectivity& c);
  void MakeVi2ti();
  void MakeLviAts();
  void MakeBaryData();
  const double* GetCoefs(size_t lvi_is_one);
  double Interp(const double c[19], const double lxi[2]) const;
  void CalcBaryTo(const double x[2], double lam[4]) const;
  void SetFnValues (size_t lvi_z1, double f_fn[3], double g_fn[3][3]) const;
  double CalcGamma(size_t vj, size_t vk) const;
  void CalcGradient(size_t lvi_1, size_t lvi_at, double grad[2]) const;
  double CalcDirectionalDeriv(const double grad[2], size_t lvi_at,
                              size_t lvi_to) const;
  void PrintState() const;
};

}}

#endif
