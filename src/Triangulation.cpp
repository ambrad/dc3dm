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
 
// It appears that our meshes are locally too uniform (e.g. b/c of RmuSmoother)
// for vertex weighting to make much of a difference.
//#define GRAD_FIT_WTS

// For debugging, make the cubic interper actually do linear.
//#define DO_LINEAR

#include <algorithm>
#include "util/include/LinAlg.hpp"
#include "Triangulation.hpp"
#include "MeshAnalyzer.hpp"

namespace util {
namespace rmesh {
using namespace std;

template<typename T> static inline ostream& operator<< (
  ostream& os, const vector<T>& v) {
  os << "[";
  for (size_t i = 0; i < v.size(); i++) {
    os << v[i];
    if (i < v.size() - 1) os << " ";
  }
  return os << "]";
}

void TriAnalyzer::Init () { InitVertexToTrisMap(); }
void TriAnalyzer::InitVertexToTrisMap () {
  _vi2ti.resize(_ma.nx());
  const size_t nt = _ma.GetTri().Size(2);
  const RectId* tri = _ma.GetTri().GetPtr();
  for (size_t ti = 1; ti <= nt; ti++, tri += 3)
    for (size_t j = 0; j < 3; j++)
      _vi2ti[tri[j]].push_back(ti);
  for (size_t i = 0; i < _vi2ti.size(); i++)
    std::sort(_vi2ti[i].begin(), _vi2ti[i].end());
}

size_t TriAnalyzer::
GetTriSharingEdge (size_t ti, RectId vi1, RectId vi2) const {
  // Use the sorted property of _vi2ti[ ].
  const vector<size_t>& tis1 = _vi2ti[vi1];
  const vector<size_t>& tis2 = _vi2ti[vi2];
  for (size_t i1 = 0, i2 = 0; i1 < tis1.size() && i2 < tis2.size(); ) {
    if (tis1[i1] == tis2[i2]) {
      if (tis1[i1] != ti) return tis1[i1];
      else { i1++; i2++; }
    } else if (tis1[i1] < tis2[i2]) i1++;
    else i2++;
  }
  return 0; // Tri ids are > 0, so 0 indicates no triangle found.
}
size_t TriAnalyzer::GetTriSharingEdge (size_t ti, size_t ei) const {
  const size_t* tri = &_ma.GetTri()(1, ti);
  return GetTriSharingEdge(ti, tri[ei], tri[(ei + 1) % 3]);
}

template<typename T>
void Uniquify (vector<T>& vi /* destroyed */, vector<T>& vo) {
  std::sort(vi.begin(), vi.end());
  typename vector<T>::iterator new_end = std::unique(vi.begin(), vi.end());
  vo.resize(new_end - vi.begin());
  for (size_t i = 0; i < vo.size(); i++) vo[i] = vi[i];
}

struct Connectivity {
  // 3 edge-sharing triangles, ordered in edge order. estis[i] is 0 if there is
  // no triangle on edge i.
  size_t estis[3];
  // The vertex-sharing tris, *including* the edge-sharing ones and ti itself,
  // in vertex order.
  vector<size_t> vstis[3];
  // Unique triangles in the support other than ti and estis.
  vector<size_t> utis;
};

CiSupport::CiSupport (const MeshAnalyzer* ma, const TriAnalyzer& ta, size_t ti)
  : _ma(*ma), _ti(ti)
{
  Connectivity c;
  SetConnectivity(ta, c);
  SetTris(c);
  MakeVi2ti();
  MakeLviAts();
  MakeBaryData();
  _cs.resize(_all_vis.size(), NULL);
  _cs_lock.resize(_cs.size());
  _cs_set.resize(_cs.size(), false);
  for (size_t i = 0; i < _cs_lock.size(); i++) omp_init_lock(&_cs_lock[i]);
}

CiSupport::~CiSupport () {
  for (size_t i = 0; i < _cs.size(); i++) {
    if (_cs[i]) delete[] _cs[i];
    omp_destroy_lock(&_cs_lock[i]);
  }
}

template<typename T> inline bool In3 (const T& e, const T set[3])
{ return e == set[0] || e == set[1] || e == set[2]; }

template<typename T> inline bool In (const T& e, const vector<T>& set) {
  return std::find(set.begin(), set.end(), e) != set.end();
}

void CiSupport::
SetConnectivity (const TriAnalyzer& ta, Connectivity& c) {
  // Make c.estis, c.vstis, _pri_vis.
  const size_t* tri = &_ma.GetTri()(1, _ti);
  _pri_vis.resize(3);
  for (size_t i = 0; i < 3; i++) {
    c.estis[i] = ta.GetTriSharingEdge(_ti, i);
    _pri_vis[i] = tri[i];
    c.vstis[i] = ta.GetTrisSharingVertex(_pri_vis[i]);
  }

  // Make _es_vis.
  for (size_t i = 0; i < 3; i++) {
    if (!c.estis[i]) continue;
    for (size_t j = 1; j <= 3; j++) {
      RectId vi = _ma.GetTri()(j, c.estis[i]);
      if (!In(vi, _pri_vis)) {
        _es_vis.push_back(vi);
        break;
      }
    }
  }

  // Make _all_vis.
  vector<RectId> av;
  for (size_t i = 0; i < 3; i++) {
    const vector<RectId>& nbrs = _ma.GetNbrs(_pri_vis[i]);
    for (size_t i = 0; i < nbrs.size(); i++) {
      RectId nvi = nbrs[i];
      if (!In(nvi, _pri_vis) && !In(nvi, _es_vis)) av.push_back(nvi);
    }
  }
  vector<RectId> uav;
  Uniquify(av, uav);
  _all_vis.reserve(3 + _es_vis.size() + uav.size());
  _all_vis = _pri_vis;
  for (size_t i = 0; i < _es_vis.size(); i++) _all_vis.push_back(_es_vis[i]);
  for (size_t i = 0; i < uav.size(); i++) _all_vis.push_back(uav[i]);

  // Make c.utis.
  vector<size_t> tis;
  tis.reserve(c.vstis[0].size() + c.vstis[1].size() + c.vstis[2].size());
  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < c.vstis[i].size(); j++) {
      size_t nti = c.vstis[i][j];
      if (nti != _ti && !In3(nti, c.estis))
        tis.push_back(nti);
    }
  Uniquify(tis, c.utis);
}

inline void Subtract2 (const double a[2], const double b[2], double c[2]) {
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
}

// Make vertex 0 of the primary triangle at (0, 0). Do all the internal
// calculations using this translation. Later, when a point to interpolate at
// comes in, translate it to the internal system.
void CiSupport::SetTris (const Connectivity& c) {
  // Map MeshAnalyzer's vertex ids to ours.
  vector<RectId> idmap(_ma.nx(), _ma.nx());
  for (size_t i = 0; i < _all_vis.size(); i++) idmap[_all_vis[i]] = i;

  // Make _t.
  { size_t n = 1;
    for (size_t i = 0; i < 3; i++) if (c.estis[i]) n++;
    _t.Resize(3, n + c.utis.size()); }
  { size_t* tri = _t.GetPtr();
    for (size_t j = 0; j < 3; j++)
      tri[j] = idmap[_pri_vis[j]];
    tri += 3;
    for (size_t i = 0; i < 3; i++) {
      if (!c.estis[i]) continue;
      for (size_t j = 0; j < 3; j++)
        tri[j] = idmap[_ma.GetTri()(j+1, c.estis[i])];
      tri += 3;
    }
    for (size_t i = 0; i < c.utis.size(); i++, tri += 3)
      for (size_t j = 0; j < 3; j++)
        tri[j] = idmap[_ma.GetTri()(j+1, c.utis[i])]; }
#ifndef NDEBUG
  for (size_t i = 1; i <= _t.Size(); i++) assert(_t(i) < _ma.nx());
#endif

  // Make _x.
  TriTag tag;
  tag.id = _ti;
  tag.dir = Dir::inside;
  tag.anchor = 0;

  _x.Resize(2, _all_vis.size());
  vector<bool> set(_all_vis.size(), false);
  // First determine the primary tri's vertices.
  { double px[3][2], origin[2];
    _ma.GetTriVertices(tag, px[0], px[1], px[2]);
    origin[0] = px[0][0]; origin[1] = px[0][1];
    for (size_t j = 0; j < 3; j++) {
      Subtract2(px[j], origin, px[j]);
      memcpy(&_x(1, j+1), px[j], 2*sizeof(double));
      set[j] = true;
    }}
  // Now use these to determine the rest.
  for (size_t i = 0; i < 3; i++) {
    double dx[2];
    Subtract2(&_x(1, i+1), &_ma.GetX()(1, _all_vis[i] + 1), dx);
    const vector<RectId>& nbrs = _ma.GetNbrs(_all_vis[i]);
    const vector<Dir::Enum>& dirs = _ma.GetNbrDirs(_all_vis[i]);
    for (size_t j = 0; j < nbrs.size(); j++) {
      assert(idmap[nbrs[j]] != _ma.nx());
      double* x = &_x(1, idmap[nbrs[j]] + 1);
      memcpy(x, &_ma.GetX()(1, nbrs[j] + 1), 2*sizeof(double));
      int Dx, Dy;
      Dir::ToVector(dirs[j], &Dx, &Dy);
      x[0] += Dx*_ma.GetDomain().dx + dx[0];
      x[1] += Dy*_ma.GetDomain().dy + dx[1];
      set[idmap[nbrs[j]]] = true;
    }
  }

  // Make the vertex adjacency graph with local vertex numbering.
  _ag.resize(_all_vis.size());
  for (size_t i = 0; i < _all_vis.size(); i++) {
    const vector<RectId>& nbrs = _ma.GetNbrs(_all_vis[i]);
    for (size_t j = 0; j < nbrs.size(); j++)
      if (idmap[nbrs[j]] != _ma.nx())
        _ag[i].push_back(idmap[nbrs[j]]);
  }

  for (size_t i = 0; i < set.size(); i++) {
    if (!set[i]) {
      cout << "i = " << i << endl;
      cout << "utis = " << c.utis << endl;
      cout << "testing: ";
      for (size_t i = 0; i < _all_vis.size(); i++)
        cout << idmap[_all_vis[i]] << " ";
      cout << endl;
      PrintState();
      assert(false);
    }
  }
}

void CiSupport::MakeVi2ti () {
  _vi2ti.resize(_all_vis.size());
  const size_t* t = _t.GetPtr();
  for (size_t ti = 1; ti <= _t.Size(2); ti++, t += 3)
    for (size_t j = 0; j < 3; j++)
      _vi2ti[t[j]].push_back(ti);
}

static inline void GetOtherLvis (size_t lvi, size_t& lvi1, size_t& lvi2) {
  lvi1 = (lvi + 1) % 3;
  lvi2 = (lvi + 2) % 3;
}

void CiSupport::MakeLviAts () {
  _lvi_ats.resize(_all_vis.size());
  for (size_t lvi = 0; lvi < _lvi_ats.size(); lvi++) {
    vector<RectId>& lvi_at = _lvi_ats[lvi];
    const vector<RectId>& nbrs = _ag[lvi];
    for (size_t i = 0; i < nbrs.size(); i++)
      if (nbrs[i] < 3) lvi_at.push_back(nbrs[i]);
    if (lvi < 3) lvi_at.push_back(lvi);
  }
}

void CiSupport::PrintState () const {
  cout << "_ti " << _ti << endl;
  cout << "_pri_vis " << _pri_vis << endl;
  cout << "_es_vis " << _es_vis << endl;
  cout << "_all_vis " << _all_vis << endl;
  cout << "_ag: " << endl;
  for (size_t i = 0; i < _ag.size(); i++)
    cout << "  " << i << " " << _ag[i] << endl;
}

void CiSupport::MakeBaryData () {
  // Primary triangle's centroid.
  const double* x1 = _x.GetPtr();
  const double* x2 = x1 + 2;
  const double* x3 = x2 + 2;
  _x4[0] = (x1[0] + x2[0] + x3[0])/3;
  _x4[1] = (x1[1] + x2[1] + x3[1])/3;

  _tri[0].Set(x2, x3, _x4);
  _tri[1].Set(x1, x3, _x4);
  _tri[2].Set(x1, x2, _x4);
  for (size_t i = 0; i < 3; i++) _tri[i].BarycentricMatrix(_Ti[i]);

#ifdef DO_LINEAR
  _tri[0].Set(x1, x2, x3);
  _tri[0].BarycentricMatrix(_Ti[0]);
#endif
}

inline bool tri_is_ccw (const double a[2], const double b[2],
                        const double c[2]) {
  // Is cross([b - a; 0], [c - a; 0]) >= 0?
  return (b[0] - a[0])*(c[1] - a[1]) >= (b[1] - a[1])*(c[0] - a[0]);
}

inline bool point_is_in_tri (const double v1[2], const double v2[2],
                             const double v3[2], const double x[2]) {
  return (tri_is_ccw(x, v1, v2) && tri_is_ccw(x, v2, v3) &&
          tri_is_ccw(x, v3, v1));
}

double CiSupport::
InterpWithOneFn (RectId vi_z1, const TriTag& tag, const double xi[2]) {
  assert(_ti == tag.id);
  // Get vertex of interest.
  int lvi_1 = -1;
  for (size_t j = 0; j < _all_vis.size(); j++)
    if (_all_vis[j] == vi_z1) {
      lvi_1 = j;
      break;
    }
  if (lvi_1 < 0) return 0;
#ifdef DO_LINEAR
  if (lvi_1 >= 3) return 0;
#endif
  // Transform xi to local coordinates.
  double lxi[2], v1[2];
  _ma.GetTriVertices(tag, v1, NULL, NULL);
  Subtract2(xi, v1, lxi);
  // Interp.
#ifdef DO_LINEAR
  return _tri[0].ToBarycentric(_Ti[0], lxi, lvi_1);
#endif
  return Interp(GetCoefs(lvi_1), lxi);
}

const size_t CiSupport::_interp_coefs[] = {
  1, 3, 3, 1, 3, 3, 1, 3, 3, 3, 3, 3, 6, 6, 6, 3, 3, 3, 1};
const size_t CiSupport::_lam_pow[][4] = {
  {3, 0, 0, 0}, {2, 1, 0, 0}, {2, 0, 1, 0}, {0, 3, 0, 0}, {1, 2, 0, 0},
  {0, 2, 1, 0}, {0, 0, 3, 0}, {1, 0, 2, 0}, {0, 1, 2, 0}, {0, 0, 2, 1},
  {0, 2, 0, 1}, {2, 0, 0, 1}, {0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1},
  {0, 0, 1, 2}, {0, 1, 0, 2}, {1, 0, 0, 2}, {0, 0, 0, 3}};

const double* CiSupport::GetCoefs (size_t lvi_z1) {
  // Quick unprotected exit. _cs_set[] is set to true only after _cs[] is given
  // the pointer, so a read of 'true' is guaranteed to give a valid pointer.
  if (_cs_set[lvi_z1]) return _cs[lvi_z1];

  omp_set_lock(&_cs_lock[lvi_z1]);
  if (!_cs[lvi_z1]) {
    // The function values (0 or 1) at the primary tri's vertices.
    double f_fn[3];
    // The direction derivatives g_fn[i][j] = d f(at vertex i)/direction(i, j)
    // at the primary tri's vertices.
    double g_fn[3][3];
    SetFnValues(lvi_z1, f_fn, g_fn);

    double* c = new double[_nc];
  
    c[c3000] = f_fn[0];
    c[c2100] = g_fn[0][1]/3 + c[c3000];
    c[c2010] = g_fn[0][2]/3 + c[c3000];
    c[c0300] = f_fn[1];
    c[c1200] = g_fn[1][0]/3 + c[c0300];
    c[c0210] = g_fn[1][2]/3 + c[c0300];
    c[c0030] = f_fn[2];
    c[c1020] = g_fn[2][0]/3 + c[c0030];
    c[c0120] = g_fn[2][1]/3 + c[c0030];
  
    c[c0021] = (c[c1020] + c[c0120] + c[c0030])/3;
    c[c0201] = (c[c1200] + c[c0300] + c[c0210])/3;
    c[c2001] = (c[c2100] + c[c2010] + c[c3000])/3;
  
    c[c0111] = 0.5*(
      CalcGamma(2, 3)*(-c[c0300] + 3*c[c0210] - 3*c[c0120] + c[c0030]) +
      (-c[c0300] + 2*c[c0210] - c[c0120] + c[c0201] + c[c0021]));
    c[c1011] = 0.5*(
      CalcGamma(3, 1)*(-c[c0030] + 3*c[c1020] - 3*c[c2010] + c[c3000]) +
      (-c[c0030] + 2*c[c1020] - c[c2010] + c[c2001] + c[c0021]));
    c[c1101] = 0.5*(
      CalcGamma(1, 2)*(-c[c3000] + 3*c[c2100] - 3*c[c1200] + c[c0300]) +
      (-c[c3000] + 2*c[c2100] - c[c1200] + c[c2001] + c[c0201]));
  
    c[c0012] = (c[c1011] + c[c0111] + c[c0021])/3;
    c[c0102] = (c[c1101] + c[c0111] + c[c0201])/3;
    c[c1002] = (c[c1101] + c[c1011] + c[c2001])/3;
    c[c0003] = (c[c1002] + c[c0102] + c[c0012])/3;

    for (size_t i = 0; i < _nc; i++) c[i] *= _interp_coefs[i];

    _cs[lvi_z1] = c;
    _cs_set[lvi_z1] = true;
  }
  omp_unset_lock(&_cs_lock[lvi_z1]);

  return _cs[lvi_z1];
}

void CiSupport::
SetFnValues (size_t lvi, double f_fn[3], double g_fn[3][3]) const {
  for (size_t i = 0; i < 3; i++) {
    f_fn[i] = 0;
    for (size_t j = 0; j < 3; j++) g_fn[i][j] = 0;
  }
  if (lvi < 3) f_fn[lvi] = 1;
  for (size_t i = 0; i < _lvi_ats[lvi].size(); i++) {
    size_t lvi1 = _lvi_ats[lvi][i], lvi2, lvi3;
    GetOtherLvis(lvi1, lvi2, lvi3);
    double grad[2];
    CalcGradient(lvi, lvi1, grad);
    g_fn[lvi1][lvi2] = CalcDirectionalDeriv(grad, lvi1, lvi2);
    g_fn[lvi1][lvi3] = CalcDirectionalDeriv(grad, lvi1, lvi3);  
  }
#ifndef NDEBUG
  for (size_t i = 0; i < 3; i++) assert(g_fn[i][i] == 0);
#endif
}

static inline double Dot2 (const double a[2], const double b[2]) {
  return a[0]*b[0] + a[1]*b[1]; }

static inline double Dot3 (const double a[3], const double b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

double CiSupport::CalcGamma (size_t vj, size_t vk) const {
  double ej4[2], ejk[2];
  const double* xj = &_x(1, vj);
  const double* xk = &_x(1, vk);
  ej4[0] = _x4[0] - xj[0];
  ej4[1] = _x4[1] - xj[1];
  ejk[0] = xk[0] - xj[0];
  ejk[1] = xk[1] - xj[1];
  return -Dot2(ej4, ejk) / Dot2(ejk, ejk);
}

inline double CiSupport::
CalcDirectionalDeriv (const double grad[2], size_t lvi_at, size_t lvi_to) const
{
  double dx[2];
  Subtract2(&_x(1, lvi_to+1), &_x(1, lvi_at+1), dx);
  return Dot2(grad, dx);
}

static inline double pow (const double lam, const size_t e) {
  switch (e) {
  case 0: return 1;
  case 1: return lam;
  case 2: return lam*lam;
  case 3: return lam*lam*lam;
  default: assert(false); return 0;
  }
}

inline double CiSupport::
Interp (const double c[19], const double lxi[2]) const {
  double lam[4];
  CalcBaryTo(lxi, lam);
  double fi = 0;
  for (size_t i = 0; i < _nc; i++) {
    double lp = c[i];
    for (size_t j = 0; j < 4; j++) lp *= pow(lam[j], _lam_pow[i][j]);
    fi += lp;
  }
  return fi;
}

inline void CiSupport::CalcBaryTo (const double x[2], double lam[4]) const {
  double lamc[3], lamk[3], mv = -1;
  size_t mti = 0;
  for (size_t i = 0; i < 3; i++) {
    _tri[i].ToBarycentric(_Ti[i], x, lamc);
    double mvc = min(min(lamc[0], lamc[1]), lamc[2]);
    // 2 of the 3 lamda vecs will have a nonpositive value. I arrange the
    // computation in this way so that *some* lambda is selected. (Otherwise,
    // finite precision could make it seem like x is not inside the tri.)
    if (i == 0 || mvc > mv) {
      mti = i;
      mv = mvc;
      memcpy(&lamk[0], &lamc[0], 3*sizeof(lamc[0]));
    }
  }
  size_t k = 0;
  for (size_t i = 0; i < mti; i++) lam[i] = lamk[k++];
  lam[mti] = 0;
  for (size_t i = mti + 1; i < 4; i++) lam[i] = lamk[k++];
}

static inline void Normalize2 (double x[2]) {
  double nrm = sqrt(x[0]*x[0] + x[1]*x[1]);
  x[0] /= nrm; x[1] /= nrm;
}

static inline void
OrganizeVertices (const size_t* t, size_t vi1, size_t& vi2, size_t& vi3) {
  size_t vi1_j = 3;
  for (size_t j = 0; j < 3; j++)
    if (t[j] == vi1) { vi1_j = j; break; }
  assert(vi1_j < 3);
  vi2 = t[(vi1_j + 1) % 3];
  vi3 = t[(vi1_j + 2) % 3];
}

static inline void SetRhs (size_t b_1, size_t col, Matrix<double>& A) {
  if (b_1 > 0) {
    for (size_t i = 1; i <= A.Size(1); i++) A(i, col) = 0;
    A(b_1, col) = 1;
  } else
    for (size_t i = 1; i <= A.Size(1); i++) A(i, col) = -1;
}

static void SetWeights (const Matrix<double>& A, Matrix<double>& oosrw) {
  // Columns [end-2 end-1] contain dx and dy.
  size_t dxc = A.Size(2) - 2, nr = A.Size(1);
  const double* x = A.GetPtr() + (dxc - 1)*nr;
  const double* y = x + nr;
  oosrw.Resize(nr);
  double* w = oosrw.GetPtr();

  double R = 0;
  for (size_t i = 0; i < nr; i++) {
    w[i] = sqrt(x[i]*x[i] + y[i]*y[i]);
    if (w[i] > R) R = w[i];
  }
  R *= 2;

  for (size_t i = 0; i < nr; i++)
    w[i] = 1/sqrt(1/w[i] - 1/R);
}

static void ApplyWeights (const Matrix<double>& oosrw, Matrix<double>& A) {
  size_t nr = A.Size(1), nc = A.Size(2);
  double* pA = A.GetPtr();
  const double* w = oosrw.GetPtr();
  for (size_t j = 0; j < nc; j++) {
    for (size_t i = 0; i < nr; i++) pA[i] *= w[i];
    pA += nr;
  }
}

void CiSupport::
CalcGradient (size_t lvi_1, size_t lvi_at, double grad[2]) const {
  assert(lvi_at < 3);
  const vector<RectId>& nbrs = _ag[lvi_at];
  const double* x_at = &_x(1, lvi_at + 1);
  if (nbrs.size() >= 5) {
    // Fit
    //     f(x,y) = z_at + c20 (x - x_at)^2 + c11 (x - x_at) (y - y_at)
    //              + c02 (y - y_at)^2 + c10 (x - x_at) + c01 (y - y_at)
    // by solving
    //     chol(W)' \ A(:,1:5) =LS chol(W)' \ A(:,1:6).
    // A's cols are ordered c20, c11, c02, c10, c01, b.

    // Form l.h.s. A(:,1:5).
    Matrix<double> A(nbrs.size(), 6);
    size_t b_1 = 0;
    for (size_t i = 0; i < nbrs.size(); i++) {
      double dx[2];
      Subtract2(&_x(1, nbrs[i]+1), x_at, dx);
      A(i+1, 1) = dx[0]*dx[0];
      A(i+1, 2) = dx[0]*dx[1];
      A(i+1, 3) = dx[1]*dx[1];
      A(i+1, 4) = dx[0];
      A(i+1, 5) = dx[1];
      if (nbrs[i] == lvi_1) b_1 = i + 1;
    }
    assert(lvi_1 == lvi_at || b_1 > 0);

    SetRhs(b_1, 6, A);

#ifdef GRAD_FIT_WTS
    if (nbrs.size() > 5) {
      Matrix<double> oosrw;
      SetWeights(A, oosrw);
      cout << oosrw << endl;
      ApplyWeights(oosrw, A);
    }
#endif

    // Factorize and solve.
    Matrix<double> R;
    SkinnyQr(A, R);
    // Calculate the gradient of f by solving
    //    R(1:5,1:5) c = R(1:5,6)
    // for just c(4:5).
    grad[1] = R(5,6) / R(5,5);
    grad[0] = (R(4,6) - R(4,5)*grad[1]) / R(4,4);
  } else {
    assert(nbrs.size() >= 2);
    // Fit just the linear part.
    Matrix<double> A(nbrs.size(), 3);
    size_t b_1 = 0;
    for (size_t i = 0; i < nbrs.size(); i++) {
      double dx[2];
      Subtract2(&_x(1, nbrs[i]+1), x_at, dx);
      A(i+1, 1) = dx[0];
      A(i+1, 2) = dx[1];
      if (nbrs[i] == lvi_1) b_1 = i + 1;
    }
    assert(lvi_1 == lvi_at || b_1 > 0);

    SetRhs(b_1, 3, A);

#ifdef GRAD_FIT_WTS
    if (_do_wts && nbrs.size() > 2) {
      Matrix<double> oosrw;
      SetWeights(A, oosrw);
      ApplyWeights(oosrw, A);
    }
#endif

    // Factorize and solve.
    Matrix<double> R;
    SkinnyQr(A, R);
    grad[1] = R(2,3) / R(2,2);
    grad[0] = (R(1,3) - R(1,2)*grad[1]) / R(1,1);
  }
}
}}
