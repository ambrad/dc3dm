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
 
#include "util/include/OpenMP.hpp"
#include "MeshAnalyzer.hpp"

namespace util {
namespace rmesh {
using namespace std;

MeshAnalyzer::MeshAnalyzer (const RectMeshUnstruct* rmu, const RmuAnalyzer* ra)
  : _rmu(rmu), _ra(ra)
{ Init(); }

void MeshAnalyzer::ResetUnitFactor (size_t uf) { Init(uf); }

void MeshAnalyzer::Init (size_t uf) {
  const vector<Rect>& rs = GetRects();
  _nr = rs.size();
  _nx = GetX().Size(2);
  // Get the unit cell size.
  _ux = rs[0].dx;
  _uy = rs[0].dy;
  for (size_t i = 1; i < rs.size(); i++) {
    _ux = std::min(_ux, rs[i].dx);
    _uy = std::min(_uy, rs[i].dy);
  }
  _ux /= uf;
  _uy /= uf;
  // Now determine the size of a cell in terms of the unit cell. In these
  // units, the cell is a square.
  _nc.resize(rs.size());
  for (size_t i = 0; i < rs.size(); i++) {
    int nx = UnitDx(rs[i].dx);
    assert (nx == UnitDy(rs[i].dy) && nx >= 1);
    _nc[i] = nx;
  }
}

void LinearInterpolatorMatrix::
Call (size_t id_z1, const vector<double>& xi, const vector<double>& yi,
      const vector<TriTag>& in, vector<double>& zi) const {
  assert (id_z1 < _ma.nx());
  for (size_t i = 0; i < xi.size(); i++) {
    const RectId* tri = &_ma.GetTri()(1, in[i].id);
    // Get vertex of interest.
    int i_tri = -1;
    for (size_t j = 0; j < 3; j++)
      if (tri[j] == id_z1) {
        i_tri = j;
        break;
      }
    if (i_tri < 0) {
      zi[i] = 0;
      continue;
    }
    // Interp. Vertex i_tri is 1 and the others are 0, so only lambda[i_tri] is
    // needed.
    double v1[2], v2[2], v3[2];
    _ma.GetTriVertices(in[i], v1, v2, v3);
    Tri2 t2(v1, v2, v3);
    double Ti[4];
    t2.BarycentricMatrix(Ti);
    double xy[2]; xy[0] = xi[i]; xy[1] = yi[i];
    zi[i] = t2.ToBarycentric(Ti, xy, i_tri);

    // Let's leave this assert output here.
#ifndef NDEBUG
    if (zi[i] <= -1.e-14) {
      const Matrix<double>& x = _ma.GetX();
      const Matrix<Dir::Enum>& tb = _ma.GetTriBdy();
      RectId rid = _ma.GetRectId(xy[0], xy[1]);
      printf("\n(%1.2f %1.2f, %1.2f %1.2f, %1.2f %1.2f) ->\n"
             "(%1.2f %1.2f, %1.2f %1.2f, %1.2f %1.2f)\n"
             "with (x, y) = (%1.2f %1.2f) and\n"
             "tri_bdy = [%d %d %d] and\n"
             "rect id = %d\n"
             "tag = id %d dir %d anchor %d\n",
             x(1, tri[0]+1), x(2, tri[0]+1),
             x(1, tri[1]+1), x(2, tri[1]+1),
             x(1, tri[2]+1), x(2, tri[2]+1),
             v1[0], v1[1], v2[0], v2[1], v3[0], v3[1],
             xy[0], xy[1],
             (int) tb(1,in[i].id), (int) tb(2,in[i].id), (int) tb(3,in[i].id),
             (int) rid,
             (int) in[i].id, (int) in[i].dir, (int) in[i].anchor);
    }
#endif
    assertpr(zi[i] >= -1.e-14, "%1.3e", zi[i]);
  }
}

void LinearInterpolatorMatrix::
GetSupportingVertices (size_t ti, vector<RectId>& vis) const {
  vis.resize(3);
  for (size_t i = 0; i < 3; i++) vis[i] = _ma.GetTri()(i+1, ti);
}

CubicInterpolatorMatrix::
CubicInterpolatorMatrix (const MeshAnalyzer* ma)
  : InterpolatorMatrix(ma), _ta(ma), _ss(ma->GetTri().Size(2), NULL)
{
  _ss_lock.resize(_ss.size());
  for (size_t i = 0; i < _ss.size(); i++)
    omp_init_lock(&_ss_lock[i]);
  _ss_set.resize(_ss.size(), false);
}
CubicInterpolatorMatrix::~CubicInterpolatorMatrix () {
  for (size_t i = 0; i < _ss.size(); i++) {
    if (_ss[i]) delete _ss[i];
    omp_destroy_lock(&_ss_lock[i]);
  }
}

void CubicInterpolatorMatrix::
Call (size_t id_z1, const vector<double>& xi, const vector<double>& yi,
      const vector<TriTag>& in, vector<double>& zi) const {
  assert (id_z1 < _ma.nx());
  for (size_t i = 0; i < xi.size(); i++) {
    InitSupportIfNeeded(in[i].id);
    double x[2];
    x[0] = xi[i]; x[1] = yi[i];
    zi[i] = _ss[in[i].id - 1]->InterpWithOneFn(id_z1, in[i], x);
  }
}

void CubicInterpolatorMatrix::
GetSupportingVertices (size_t ti, vector<RectId>& vis) const {
  InitSupportIfNeeded(ti);
  _ss[ti-1]->GetSupportingVertices(ti, vis);
}

inline void CubicInterpolatorMatrix::InitSupportIfNeeded (size_t ti) const {
  // Quick unprotected exit.
  if (_ss_set[ti - 1]) return;
  // Deal with parallel access if we need a new CiSupport.
  omp_set_lock(&_ss_lock[ti - 1]);
  if (!_ss[ti - 1]) {
    CiSupport* cis = new CiSupport(&_ma, _ta, ti);
    _ss[ti - 1] = cis;
    _ss_set[ti - 1] = true;
  }
  omp_unset_lock(&_ss_lock[ti - 1]);
}

}}
