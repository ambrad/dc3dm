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
#include <algorithm>
#include <queue>
#include <set>
#include <list>
#include <limits>
#include "util/include/Util.hpp"
#include "util/include/CodeAnalysis.hpp"
#include "util/include/Debug.hpp"
#include "util/include/LinAlg.hpp"
#include "Tri2.hpp"
#include "RectMeshUnstruct_pri.hpp"

namespace util {
namespace rmesh {
static const size_t nkids = 4;

// -----------------------------------------------------------------------------
// Plane geometry.

template<typename T>
struct PointT {
  T x, y;
  PointT () : x(0), y(0) {}
  PointT (T ix, T iy) : x(ix), y(iy) {}
  void Set (T ix, T iy) { x = ix; y = iy; }
};

template<typename T>
class PointLineSegT {
public:
  PointLineSegT () {}
  PointLineSegT(const PointT<T>& ip, const PointT<T>& iq) : _p(ip), _q(iq) {}
  PointT<T>& p () { return _p; }
  PointT<T>& q () { return _q; }
  const PointT<T>& p () const { return _p; }
  const PointT<T>& q () const { return _q; }
  static bool IntersectOpen(const PointLineSegT<T>& ls1,
                            const PointLineSegT<T>& ls2);
private:
  PointT<T> _p, _q;
};

template<typename T>
class PointRectT {
public:
  PointRectT () {}
  PointRectT (const PointT<T>& c1, const PointT<T>& c2)
    : _ll(c1), _ur(c2) { SetCorners(); }
  PointRectT (const PointLineSegT<T>& ls) // Rect spanning the PLS.
    : _ll(ls.p()), _ur(ls.q()) { SetCorners(); }
  PointRectT (const RectT<T>& r)
    : _ll(r.x, r.y), _ur (r.x + r.dx, r.y + r.dy) {}
  PointRectT (const RectT<T>& r, T dx, T dy)
    : _ll(r.x + dx, r.y + dy), _ur (r.x + r.dx + dx, r.y + r.dy + dy) {}
  PointT<T>& ll () { return _ll; } // Lower-left corner.
  PointT<T>& ur () { return _ur; } // Upper-right corner.
  const PointT<T>& ll () const { return _ll; }
  const PointT<T>& ur () const { return _ur; }

  // Intersection of open rectangles.
  static bool IntersectOpen (const PointRectT<T>& r1, const PointRectT<T>& r2) {
    return !(r1.ur().x <= r2.ll().x || r1.ll().x >= r2.ur().x ||
             r1.ur().y <= r2.ll().y || r1.ll().y >= r2.ur().y);
  }
  // Intersection of closed rectangles.
  static bool IntersectClosed (const PointRectT<T>& r1,
                               const PointRectT<T>& r2) {
    return !(r1.ur().x < r2.ll().x || r1.ll().x > r2.ur().x ||
             r1.ur().y < r2.ll().y || r1.ll().y > r2.ur().y);
  }

  static void Distance(const PointRectT<T>& r1, const PointRectT<T>& r2,
                       T* dx, T* dy);
  static void Distance(const PointRectT<T>& r, const PointT<T>& p,
                       T* dx, T* dy);
private:
  PointT<T> _ll, _ur;

  void SetCorners ()
  { if (_ll.x > _ur.x) std::swap(_ll.x, _ur.x);
    if (_ll.y > _ur.y) std::swap(_ll.y, _ur.y); }
};

// Preferred types.
typedef PointT<double> Point;
typedef PointLineSegT<int> PointLineSeg;
typedef PointRectT<double> PointRect;

template<typename T>
ostream& operator<< (ostream& os, const PointLineSegT<T>& pls) {
  return os << "[(" << pls.p().x << ", " << pls.p().y << ") to (" << pls.q().x
            << ", " << pls.q().y << ")]";
}

template<typename T>
ostream& operator<< (ostream& os, const RectT<T>& r) {
  return os << "[(" << r.x << ", " << r.y << ") ("
            << r.x + r.dx << ", " << r.y + r.dy << ")]";
}

template<typename T>
ostream& operator<< (ostream& os, const PointRectT<T>& r) {
  return os << "[(" << r.ll().x << ", " << r.ll().y << ") ("
            << r.ur().x << ", " << r.ur().y << ")]";
}

template<typename T>
inline bool operator== (const PointT<T>& p1, const PointT<T>& p2)
{ return p1.x == p2.x && p1.y == p2.y; }

template<typename T>
inline void PointRectT<T>::
Distance (const PointRectT<T>& r1, const PointRectT<T>& r2, T* dx, T* dy) {
  if (dx)
    *dx = (r1.ll().x <= r2.ll().x) ?
      std::max((T) 0, r2.ll().x - r1.ur().x) :
      std::max((T) 0, r1.ll().x - r2.ur().x);
  if (dy)
    *dy = (r1.ll().y <= r2.ll().y) ?
      std::max((T) 0, r2.ll().y - r1.ur().y) :
      std::max((T) 0, r1.ll().y - r2.ur().y);
}

template<typename T>
inline void PointRectT<T>::
Distance (const PointRectT<T>& r, const PointT<T>& p, T* dx, T* dy) {
  if (dx)
    *dx = (r.ll().x <= p.x) ?
      std::max((T) 0, p.x - r.ur().x) :
      std::max((T) 0, r.ll().x - p.x);
  if (dy)
    *dy = (r.ll().y <= p.y) ?
      std::max((T) 0, p.y - r.ur().y) :
      std::max((T) 0, r.ll().y - p.y);
}

template<typename T>
void RectT<T>::Distance (const RectT<T>& r1, const RectT<T>& r2, T* dx, T* dy) {
  PointRectT<T>::Distance(PointRectT<T>(r1), PointRectT<T>(r2), dx, dy);
}

// Intersection of open line segments.
template<typename T>
bool PointLineSegT<T>::
IntersectOpen (const PointLineSegT<T>& ls1, const PointLineSegT<T>& ls2) {
  // Quick exit on rect intersection test.
  if (!PointRectT<T>::IntersectOpen(PointRectT<T>(ls1), PointRectT<T>(ls2)))
    return false;
  // Quick exit on line segments being ordered oppositely.
  if (ls1.p() == ls2.q() && ls1.q() == ls2.p()) return true;
  // (We could also test for the same orientation, but in our application,
  // only oppositely ordered happens.)
  // if (ls1.q() == ls2.q() && ls1.p() == ls2.p()) return true;
  // Quick exit on line segments meeting at their endpoints. We already know
  // *both* end points do not.
  if (ls1.p() == ls2.p() || ls1.p() == ls2.q() ||
      ls1.q() == ls2.p() || ls1.q() == ls2.q()) return false;
  // Form the equation corresponding to the intersection of parametric lines:
  //     ls1.p + (ls1.q - ls1.p) w(1) = ls2.p + (ls2.q - ls2.p) w(2).
  // This yields A w = b, and we solve for w. If both w are in [0,1], then the
  // segments intersect.
  double a11, a12, a21, a22, b1, b2, det, w1, w2;
  a11 = ls1.q().x - ls1.p().x; a21 = ls1.q().y - ls1.p().y;
  a12 = ls2.p().x - ls2.q().x; a22 = ls2.p().y - ls2.q().y;
  det = a11*a22 - a12*a21;
  if (det == 0) return false;
  b1 = ls2.p().x - ls1.p().x; b2 = ls2.p().y - ls1.p().y;
  w1 = ( a22*b1 - a12*b2) / det;
  if (w1 <= 0 || w1 >= 1) return false;
  w2 = (-a21*b1 + a11*b2) / det;
  return w2 > 0 && w2 < 1;
}

// -----------------------------------------------------------------------------
// Obtain an f-satisfying mesh.

const double RectOpts::max_aspect_ratio = 1.1;
RectOpts::RectOpts (FILE* fid) { Deserialize(fid); }
void RectOpts::Serialize (FILE* fid) const {
  write(&min_len, 2, fid);
  write(&uniform, 1, fid);
}
void RectOpts::Deserialize (FILE* fid) {
  read(&min_len, 2, fid);
  read(&uniform, 1, fid);
}

void DeserializeRect (Rect& r, FILE* fid) { read(&r.x, 4, fid); }
void SerializeRect (const Rect& r, FILE* fid) { write(&r.x, 4, fid); }
void InitRect (Rect& r, FILE* fid) { DeserializeRect(r, fid); }

bool operator< (const Rect& r1, const Rect& r2) {
  if (r1.dx < r2.dx) return true;
  if (r1.dx == r2.dx) return r1.dy < r2.dy;
  return false;
}
    
// i is a base-1 index.
static inline void
GetSurroundingPoints (const Matrix<double>& v, double x, size_t* i) {
  const double* pl = std::lower_bound(v.GetPtr(), v.GetPtr() + v.Size(), x);
  *i = (size_t) (pl - v.GetPtr());
}

void TensorMeshLinInterpRF::
Call (const vector<double>& x, const vector<double>& y, vector<double>& f) {
  for (size_t i = 0, n = x.size(); i < n; i++) {
    double alpha, f1, f2;
    size_t ix, iy;
    GetSurroundingPoints(_x, x[i], &ix);
    GetSurroundingPoints(_y, y[i], &iy);
    alpha = (y[i] - _y(iy)) / (_y(iy + 1) - _y(iy));
    f1 = (1 - alpha) * _f(iy, ix    ) + alpha * _f(iy + 1, ix    );
    f2 = (1 - alpha) * _f(iy, ix + 1) + alpha * _f(iy + 1, ix + 1);
    alpha = (x[i] - _x(ix)) / (_x(ix + 1) - _x(ix));
    f[i] = (1 - alpha) * f1 + alpha * f2;
  }
}

TensorMeshLinInterpRF::
TensorMeshLinInterpRF (const Matrix<double>& x, const Matrix<double>& y,
                       const Matrix<double>& f)
  : _x(x), _y(y), _f(f)
{
  assert(_y.Size() == _f.Size(1));
  assert(_x.Size() == _f.Size(2));
  assert(_x.Size() >= 2 && _y.Size() >= 2);
}

const Matrix<double>& TensorMeshLinInterpRF::GetX () const { return _x; }
const Matrix<double>& TensorMeshLinInterpRF::GetY () const { return _y; }
const Matrix<double>& TensorMeshLinInterpRF::GetF () const { return _f; }

namespace getmin {
inline size_t Low (const Matrix<double>& x, double v) {
  const double* xlo =
    std::lower_bound(x.GetPtr(), x.GetPtr() + x.Size(), v);
  if (*xlo > v && xlo > x.GetPtr()) xlo--;
  return xlo - x.GetPtr() + 1;
}

inline size_t High (const Matrix<double>& x, size_t xlo, double v) {
  return std::lower_bound(x.GetPtr() + xlo, x.GetPtr() + x.Size(), v)
    - x.GetPtr() + 1;
}
}

double TensorMeshLinInterpRF::GetMin (const Rect& dmn) const {
  using namespace getmin;
  size_t jlo = Low(_x, dmn.x);
  // min is to prevent numerical error from setting jhi to _x.Size() + 1.
  size_t jhi = std::min(_x.Size(), High(_x, jlo, dmn.x + dmn.dx));
  size_t ilo = Low(_y, dmn.y);
  size_t ihi = std::min(_y.Size(), High(_y, ilo, dmn.y + dmn.dy));
  assertpr(jlo >= 1 && jhi <= _f.Size(2) && ilo >= 1 && ihi <= _f.Size(1),
           "jlo %d jhi %d ilo %d ihi %d f.Size() = %d %d", (int) jlo, (int) jhi,
           (int) ilo, (int) ihi, (int) _f.Size(1), (int) _f.Size(2));
  double f_min = _f(ilo, jlo);
  for (size_t j = jlo; j <= jhi; j++)
    for (size_t i = ilo; i <= ihi; i++)
      f_min = std::min(f_min, _f(i,j));
  return f_min;
}

void TensorMeshLinInterpRF::GetMin (double* x, double* y, double* f) const {
  *f = _f(1, 1), *x = _x(1), *y = _y(1);
  for (size_t j = 1; j <= _y.Size(); j++)
    for (size_t i = 1; i <= _x.Size(); i++) {
      if (_f(i,j) < *f) {
        *f = _f(i,j);
        *x = _x(j);
        *y = _y(i);
      }
    }
}

RectMeshUnstruct::QuadTree::QuadTree (const Rect& r)
  : _r(r), _is_leaf(true) {}
    
RectMeshUnstruct::QuadTree::QuadTree (FILE* fid)
  : _is_leaf(true)
{ Deserialize(fid); }

RectMeshUnstruct::QuadTree::QuadTree (const RectMeshUnstruct::QuadTree& qt)
  : _r(qt._r), _is_leaf(qt._is_leaf)
{
  if (_is_leaf)
    _l.id = qt._l.id;
  else
    for (size_t i = 0; i < nkids; i++)
      _n.kids[i] = new QuadTree(*qt._n.kids[i]);
}

RectMeshUnstruct::QuadTree::~QuadTree () {
  if (!_is_leaf)
    for (size_t i = 0; i < nkids; i++)
      if (_n.kids[i]) delete _n.kids[i];
}

void RectMeshUnstruct::QuadTree::Serialize (FILE* fid) {
  SerializeRect(_r, fid);
  char is_leaf = _is_leaf;
  write(&is_leaf, 1, fid);
  if (!_is_leaf)
    for (size_t i = 0; i < nkids; i++) _n.kids[i]->Serialize(fid);
  else
    write(&_l.id, 1, fid);
}

void RectMeshUnstruct::QuadTree::Deserialize (FILE* fid) {
  DeserializeRect(_r, fid);
  char is_leaf;
  read(&is_leaf, 1, fid);
  _is_leaf = is_leaf;
  if (!_is_leaf) {
    size_t i;
    try {
      for (i = 0; i < nkids; i++) _n.kids[i] = new QuadTree(fid);
    } catch (const FileException& e) {
      for (size_t j = 0; j < i; j++) delete _n.kids[i];
      throw e;
    }
  } else
    read(&_l.id, 1, fid);
}

void RectMeshUnstruct::QuadTree::Split () {
  double dxh = 0.5*_r.dx, dyh = 0.5*_r.dy;
  _n.kids[0] = new QuadTree(Rect(_r.x + dxh, _r.y + dyh, dxh, dyh));
  _n.kids[1] = new QuadTree(Rect(_r.x      , _r.y + dyh, dxh, dyh));
  _n.kids[2] = new QuadTree(Rect(_r.x      , _r.y      , dxh, dyh));
  _n.kids[3] = new QuadTree(Rect(_r.x + dxh, _r.y      , dxh, dyh));
  _is_leaf = false;
}

void RectMeshUnstruct::QuadTree::
Split (const RectOpts& ro, ResolutionFn* rf, RectId* id, bool first) {
  if (first ||
      max(_r.dx, _r.dy) > min(ro.max_len, max(rf->GetMin(_r), ro.min_len))) {
    Split();
    for (size_t i = 0; i < nkids; i++) _n.kids[i]->Split(ro, rf, id);
  } else {
    _l.id = *id;
    (*id)++;
  }
}

void RectMeshUnstruct::QuadTree::PushBackRects (vector<Rect>* rs) {
  if (!_is_leaf)
    for (size_t i = 0; i < nkids; i++) _n.kids[i]->PushBackRects(rs);
  else
    rs->push_back(_r);
}

inline int RectMeshUnstruct::QuadTree::GetRectId (double x, double y) const {
  // RectMeshUnstruct::GetRectId assertion protects us: (x, y) is guaranteed
  // to be in _r.
  if (!_is_leaf) {
    if (x >= _r.x + 0.5*_r.dx) {
      if (y >= _r.y + 0.5*_r.dy) return _n.kids[0]->GetRectId(x, y);
      else return _n.kids[3]->GetRectId(x, y);
    } else {
      if (y >= _r.y + 0.5*_r.dy) return _n.kids[1]->GetRectId(x, y);
      else return _n.kids[2]->GetRectId(x, y);
    }
  } else return _l.id;
}

void RectMeshUnstruct::QuadTree::Refine (RectId* id) {
  if (_is_leaf) {
    Split();
    for (size_t i = 0; i < nkids; i++) _n.kids[i]->_l.id = (*id)++;
  } else {
    for (size_t i = 0; i < nkids; i++) _n.kids[i]->Refine(id);
  }
}

static void
BreakUpDomain (const RectOpts& ro, double Dx, double Dy,
               double* dx, double* dy, size_t* nx, size_t* ny) {
  if (Dx < Dy) {
    BreakUpDomain(ro, Dy, Dx, dy, dx, ny, nx);
    return;
  }
  if (ro.uniform) {
    // The factor of 4 is because of the comment below.
    *ny = (size_t) ceil(Dy / (4 * ro.max_len));
  } else
    *ny = 1;
  for ( ; ; (*ny)++) {
    *dy = Dy / *ny;
    *nx = (size_t) round(Dx / *dy);
    *dx = Dx / *nx;
    if (*dy / *dx <= ro.max_aspect_ratio &&
        *dx / *dy <= ro.max_aspect_ratio) break;
  }
  if (*nx == 1 || *ny == 1) {
    // This, combined with first=true in split, assures that the domain has a
    // minimum of 4 rectangles in each direction. If there are fewer, in
    // particular 2 or 3 in a periodic direction, the triangle-making procedure
    // can mess up. Of course, we expect domains to be way more discretized than
    // this, so this is for safety.
    *nx *= 2; *ny *= 2;
    *dx /= 2; *dy /= 2; 
  }
}

RectMeshUnstruct::
RectMeshUnstruct (const Rect& domain, const RectOpts& ro, ResolutionFn* rf)
  : _domain(domain), _ro(ro)
{
  // First break the domain into squares.
  double dx, dy;
  BreakUpDomain(ro, _domain.dx, _domain.dy, &dx, &dy, &_nx, &_ny);
#if 0
  printf("%% domain %1.2f %1.2f nx ny %d %d dx dy %1.2f %1.2f ratio %1.2f\n",
         _domain.dx, _domain.dy, (int) _nx, (int) _ny, dx, dy,
         std::max(dx/dy, dy/dx));
#endif

  // Now split each square.
  _qts.reserve(_ny * _nx);
  RectId id = 0;
  for (size_t iy = 0, k = 0; iy < _ny; iy++)
    for (size_t ix = 0; ix < _nx; ix++, k++) {
      QuadTree* qt = new QuadTree
        (Rect(_domain.x + ix*dx, _domain.y + iy*dy, dx, dy));
      qt->Split(ro, rf, &id, true);
      _qts.push_back(qt);
    }
}

RectMeshUnstruct::
RectMeshUnstruct (const RectMeshUnstruct& rmu)
  : _domain(rmu._domain), _ro(rmu._ro), _nx(rmu._nx), _ny(rmu._ny)
{
  _qts.resize(rmu._qts.size());
  for (size_t i = 0; i < _qts.size(); i++)
    _qts[i] = new QuadTree(*rmu._qts[i]);
}

RectMeshUnstruct::RectMeshUnstruct (const string& filename)
  throw (FileException)
{
  FILE* fid = fopen(filename.c_str(), "r");
  if (!fid) throw FileException("Can't read " + filename);
  Deserialize(fid);
  fclose(fid);
}

RectMeshUnstruct::~RectMeshUnstruct ()
{ for (size_t i = 0; i < _qts.size(); i++) delete _qts[i]; }

RectMeshUnstruct* NewRectMeshUnstruct (const Rect& domain, const RectOpts& ro,
                                       ResolutionFn* rf)
{ return new RectMeshUnstruct(domain, ro, rf); }
RectMeshUnstruct* NewRectMeshUnstruct (const string& filename)
  throw (FileException)
{ return new RectMeshUnstruct(filename); }
void DeleteRectMeshUnstruct (RectMeshUnstruct* rmu) { delete rmu; }

void RectMeshUnstruct::Serialize (const string& filename) const
  throw (FileException)
{
  FILE* fid = fopen(filename.c_str(), "w");
  if (!fid) throw FileException("Can't read " + filename);
  Serialize(fid);
  fclose(fid);
}
    
void RectMeshUnstruct::Serialize (FILE* fid) const throw (FileException) {
  write("rmu ", 4, fid);
  SerializeRect(_domain, fid);
  _ro.Serialize(fid);
  write(&_nx, 1, fid);
  write(&_ny, 1, fid);
  long long int n = _qts.size();
  write(&n, 1, fid);
  for (size_t i = 0; i < _qts.size(); i++) _qts[i]->Serialize(fid);
}

void RectMeshUnstruct::Deserialize (FILE* fid) {
  try {
    char tag[4];
    read(tag, 4, fid);
    if (string(tag, 4) != string("rmu "))
      throw FileException("Not a RectMeshUnstruct serial file.");
    DeserializeRect(_domain, fid);
    _ro.Deserialize(fid);
    read(&_nx, 1, fid);
    read(&_ny, 1, fid);
    long long int n;
    read(&n, 1, fid);
    _qts.reserve(n);
    for (size_t i = 0; i < (size_t) n; i++) {
      QuadTree* qt = new QuadTree(fid);
      _qts.push_back(qt);
    }
  } catch (const FileException& e) {
    for (size_t i = 0; i < _qts.size(); i++) delete _qts[i];
    throw e;
  }
}

const Rect& RectMeshUnstruct::GetDomain () const { return _domain; }

const vector<Rect>& RectMeshUnstruct::GetRects () const {
  if (_rs.empty()) // Lazy evaluation
    for (size_t i = 0; i < _qts.size(); i++)
      _qts[i]->PushBackRects(&_rs);
  return _rs;
}

template<typename RectType, typename T>
inline void PositionToIndices 
(const RectType& r, const T x, const T y, const size_t nx,
 const size_t ny, size_t* ix, size_t* iy) {
  *ix = (size_t) ((x - r.x) / r.dx * nx);
  if (*ix == nx) (*ix)--;
  *iy = (size_t) ((y - r.y) / r.dy * ny);
  if (*iy == ny) (*iy)--;
}

int RectMeshUnstruct::GetRectId (double x, double y) const {
  if (x < _domain.x || x > _domain.x + _domain.dx ||
      y < _domain.y || y > _domain.y + _domain.dy || _qts.empty())
    return -1;
  // Get first-level rectangle and its quadtree.
  size_t ix, iy;
  PositionToIndices(_domain, x, y, _nx, _ny, &ix, &iy);
  const QuadTree* qt = _qts[_nx * iy + ix];
  // Now find the rectangle within the quadtree.
  return qt->GetRectId(x, y);
}

void RectMeshUnstruct::Refine () {
  RectId id = 0;
  for (size_t i = 0; i < _qts.size(); i++) _qts[i]->Refine(&id);
  _rs.clear();
}

// -----------------------------------------------------------------------------
// Analyze a mesh.

namespace {
  void WriteTris (const RmuAnalyzer& ra) {
    printf("Writing tris\n");
    const size_t nt = ra.GetTri().Size(2);
    Matrix<double> A(6, nt);
    double* p = A.GetPtr();
    for (size_t i = 1; i <= nt; i++) {
      TriTag tag;
      tag.id = i;
      tag.anchor = 0;
      tag.dir = Dir::inside;
      ra.GetPeriodicTri2(tag, p, p+2, p+4);
      p += 6;
    }
    util::WriteMatrix(A, "/scratch/ambrad/testrmesh/tris.mat");
  }
}

inline void GetUnitCellSize 
(const RectMeshUnstruct& rmu, double* dx, double* dy)
{
  const vector<Rect>& rs = rmu.GetRects();
  *dx = rs[0].dx;
  *dy = rs[0].dy;
  for (size_t i = 1; i < rs.size(); i++) {
    *dx = min(*dx, rs[i].dx);
    *dy = min(*dy, rs[i].dy);
  }
}

RmuAnalyzer::RmuAnalyzer (const RectMeshUnstruct* rmu, const Boundaries& b)
  : _rmu(*rmu), _rs(rmu->GetRects()), _r(rmu->GetDomain()), _b(b),
    _nr(_rs.size())
{
  GetUnitCellSize(_rmu, &_ux, &_uy);
  _any_periodic = _b.GetBC(Dir::E) == Boundaries::bc_periodic
    || _b.GetBC(Dir::N) == Boundaries::bc_periodic;

  Timer* t = Ca::GetTimer();
  double et_ois, et_t;
  { Vve age;
    t->Tic(1);
    OrganizeInSpace(age);
    et_ois = t->Toc(1); t->Tic(1);

    // Triangulation.
    Matrix<Edge> trie;
    t->Tic(1);
    Triangulate(age, trie);
    et_t = t->Toc(1); t->Tic(1);

    // Convert to final data structures.
    _ag.resize(age.size());
    _ag_bdy.resize(age.size());
    for (size_t i = 0; i < age.size(); i++) {
      _ag[i].resize(age[i].size());
      _ag_bdy[i].resize(age[i].size());
      for (size_t j = 0; j < age[i].size(); j++) {
        _ag[i][j] = age[i][j].id;
        _ag_bdy[i][j] = age[i][j].bdy;
      }
    }
    _tri.Resize(3, trie.Size(2));
    _tri_bdy.Resize(3, trie.Size(2));
    for (size_t j = 1; j <= trie.Size(2); j++)
      for (size_t i = 1; i <= 3; i++) {
        _tri(i,j) = trie(i,j).id;
        _tri_bdy(i,j) = trie(i,j).bdy;
      }
    //WriteTris(*this);
  }
  double et_c = t->Toc(1); t->Tic(1);

  // For fast triangle finding.
  MakeRtoT();
  double et_r2t = t->Toc(1);
  caprint("RA: ois (%1.2e) 1.00 tri %1.2f convert %1.2f r2t %1.2f\n",
          et_ois, et_t/et_ois, et_c/et_ois, et_r2t/et_ois);

#ifndef NDEBUG
  // Assertions about the symmetry proprties of _ag and _ag_bdy, as well as
  // uniqueness and sortedness of edge lists.
  { size_t nempty = 0;
    for (size_t i = 0; i < _ag.size(); i++) {
      if (_ag[i].empty()) {
        nempty++;
        assert(Dir::IsCorner(GetBoundaryDir(i)));
        continue;
      }
      for (size_t j = 0; j < _ag[i].size(); j++) {
        // Unique and sorted.
        if (j > 0) assert(_ag[i][j] > _ag[i][j-1]);
        // Symmetry.
        RectId agij = _ag[i][j];
        bool fnd = false;
        for (size_t k = 0; k < _ag[agij].size(); k++)
          if (_ag[agij][k] == i) {
            // (1) In direction pointing.
            assert(_ag_bdy[agij][k] == Dir::Opposite(_ag_bdy[i][j]));
            fnd = true;
            break; }
        // (2) Every edge has an associated oppposite edge.
        assert(fnd); }}
    assert(nempty == 0 || nempty == 2); }
#endif
}

RmuAnalyzer::RmuAnalyzer (const RectMeshUnstruct* rmu, const Boundaries& b,
                          bool want_ag_only)
  : _rmu(*rmu), _rs(rmu->GetRects()), _r(rmu->GetDomain()), _b(b),
    _nr(_rs.size())
{
  assert(want_ag_only);
  GetUnitCellSize(_rmu, &_ux, &_uy);
  _any_periodic = _b.GetBC(Dir::E) == Boundaries::bc_periodic
    || _b.GetBC(Dir::N) == Boundaries::bc_periodic;

  Vve age;
  OrganizeInSpace(age);

  _ag.resize(age.size());
  _ag_bdy.resize(age.size());
  for (size_t i = 0; i < age.size(); i++) {
    _ag[i].resize(age[i].size());
    _ag_bdy[i].resize(age[i].size());
    for (size_t j = 0; j < age[i].size(); j++) {
      _ag[i][j] = age[i][j].id;
      _ag_bdy[i][j] = age[i][j].bdy;
    }
  }
}

RmuAnalyzer::RmuAnalyzer (const RectMeshUnstruct* rmu, const string& filename)
  throw (FileException)
  : _rmu(*rmu), _rs(rmu->GetRects()), _r(rmu->GetDomain())
{
  FILE* fid = fopen(filename.c_str(), "r");
  if (!fid) throw FileException("Can't read " + filename);
  Deserialize(fid);
  fclose(fid);
  if (_nr != _rs.size()) throw FileException("rmu->GetRects().size() != _nr");
}

RmuAnalyzer::~RmuAnalyzer () {}

void RmuAnalyzer::Serialize (const string& filename) const throw (FileException)
{
  FILE* fid = fopen(filename.c_str(), "w");
  if (!fid) throw FileException("Can't read " + filename);
  Serialize(fid);
  fclose(fid);
}

void RmuAnalyzer::Serialize (FILE* fid) const throw (FileException) {
  size_t n;

  write("ra  ", 4, fid);
  _b.Serialize(fid);
  write(&_ux, 1, fid);
  write(&_uy, 1, fid);
  write(&_nr, 1, fid);
  write(&_any_periodic, 1, fid);

  n = _ag.size(); write(&n, 1, fid);
  for (size_t i = 0; i < n; i++) Write(_ag[i], fid);

  n = _ag_bdy.size();
  write(&n, 1, fid);
  for (size_t i = 0; i < n; i++) Write(_ag_bdy[i], fid);

  _x.Serialize(fid);  
  Write(_bdypt_side, fid);
  _tri.Serialize(fid);
  _tri_bdy.Serialize(fid);

  n = _r2t.size();
  write(&n, 1, fid);
  for (size_t i = 0; i < n; i++) Write(_r2t[i], fid);
}

void RmuAnalyzer::Deserialize (FILE* fid) {
  size_t n;

  char tag[4];
  read(tag, 4, fid);
  if (string(tag, 4) != string("ra  ")) {
    printf("tag = %s\n", string(tag, 4).c_str());
    throw FileException("Not an RmuAnalyzer serial file.");
  }

  _b.Deserialize(fid);
  read(&_ux, 1, fid);
  read(&_uy, 1, fid);
  read(&_nr, 1, fid);
  read(&_any_periodic, 1, fid);

  read(&n, 1, fid); _ag.resize(n);
  for (size_t i = 0; i < n; i++) Read(_ag[i], fid);

  read(&n, 1, fid); _ag_bdy.resize(n);
  for (size_t i = 0; i < n; i++) Read(_ag_bdy[i], fid);

  _x.Deserialize(fid);
  Read(_bdypt_side, fid);
  _tri.Deserialize(fid);
  _tri_bdy.Deserialize(fid);

  read(&n, 1, fid); _r2t.resize(n);
  for (size_t i = 0; i < n; i++) Read(_r2t[i], fid);
}

RmuAnalyzer* NewRmuAnalyzer (const RectMeshUnstruct* rmu, const Boundaries& b)
{ return new RmuAnalyzer(rmu, b); }
RmuAnalyzer* NewRmuAnalyzer (const RectMeshUnstruct* rmu,
                             const string& filename)
  throw (FileException)
{ return new RmuAnalyzer(rmu, filename); }
void DeleteRmuAnalyzer (RmuAnalyzer* ra) { delete ra; }

// -----------------------------------
// Adjacency graph.

Dir::Enum RmuAnalyzer::OnWhichBoundary (const Rect& r) const {
  const double hux = 0.5*_ux, huy = 0.5*_uy;
  if        (r.x        - hux < _r.x        ) {
    if      (r.y        - huy < _r.y        ) return Dir::SW;
    else if (r.y + r.dy + huy > _r.y + _r.dy) return Dir::NW;
    else                                      return Dir::W;
  } else if (r.x + r.dx + hux > _r.x + _r.dx) {
    if      (r.y        - huy < _r.y        ) return Dir::SE;
    else if (r.y + r.dy + huy > _r.y + _r.dy) return Dir::NE;
    else                                      return Dir::E;
  } else if (r.y        - huy < _r.y        ) return Dir::S;
  else if   (r.y + r.dy + huy > _r.y + _r.dy) return Dir::N;
  else                                        return Dir::inside;
}

Dir::Enum RmuAnalyzer::OnWhichBoundary (const double x, const double y) const {
  const double hux = 0.5*_ux, huy = 0.5*_uy;
  if        (x - hux < _r.x        ) {
    if      (y - huy < _r.y        ) return Dir::SW;
    if      (y + huy > _r.y + _r.dy) return Dir::NW;
    else                             return Dir::W;
  } else if (x + hux > _r.x + _r.dx) {
    if      (y - huy < _r.y        ) return Dir::SE;
    if      (y + huy > _r.y + _r.dy) return Dir::NE;
    else                             return Dir::E;
  } else if (y - huy < _r.y        ) return Dir::S;
  else if   (y + huy > _r.y + _r.dy) return Dir::N;
  else                               return Dir::inside;
}

// The smallest size in unit coords is 2, not 1, so that the middle of the cell
// in unit coords is an integer.
inline int RmuAnalyzer::UnitDx (double dx) const
{ return (int) round(2 * dx / _ux); }
inline int RmuAnalyzer::UnitDy (double dy) const
{ return (int) round(2 * dy / _uy); }
inline int RmuAnalyzer::UnitX (double x) const { return UnitDx(x - _r.x); }
inline int RmuAnalyzer::UnitY (double y) const { return UnitDy(y - _r.y); }
inline double RmuAnalyzer::DomainX (int x) const
{ return _r.x + 0.5 * _ux * x; }
inline double RmuAnalyzer::DomainY (int y) const
{ return _r.y + 0.5 * _uy * y; }

inline bool InsideRect (const Rect& r, double x, double y)
{ return x >= r.x && x < r.x + r.dx && y >= r.y && y < r.y + r.dy; }

inline bool InsideRect (const Rect& r_out, const Rect& r_in) {
  return r_in.x >= r_out.x && r_in.x + r_in.dx <= r_out.x + r_out.dx &&
    r_in.y >= r_out.y && r_in.y + r_in.dy <= r_out.y + r_out.dy;
}

inline bool NoPeriodic (const Boundaries& b) {
  return b.GetBC(Dir::E) != Boundaries::bc_periodic
    &&   b.GetBC(Dir::N) != Boundaries::bc_periodic;
}

// Set (x, y) on output to the periodic coordinate in GetDomain(). Return false
// if not possible. Unchanged if (x, y) in GetDomain() on input. Set bdy to the
// boundary one goes through get from the input (x, y) to the output (x, y);
// equivalently, to the direction one goes to get from the output to the input.
bool RmuAnalyzer::MapToDomain (double* x, double* y, Dir::Enum* bdy) const {
  const Rect& d = _rmu.GetDomain();
  *bdy = Dir::inside;
  if (InsideRect(d, *x, *y)) return true;
  if (*x > d.x + d.dx) {
    if (_b.GetBC(Dir::E) != Boundaries::bc_periodic) return false;
    *bdy = Dir::E;
    *x -= d.dx;
  } else if (*x < d.x) {
    if (_b.GetBC(Dir::W) != Boundaries::bc_periodic) return false;
    *bdy = Dir::W;
    *x += d.dx;
  }
  if (*y > d.y + d.dy) {
    if (_b.GetBC(Dir::N) != Boundaries::bc_periodic) return false;
    switch (*bdy) {
    case Dir::inside: *bdy = Dir::N; break;
    case Dir::E: *bdy = Dir::NE; break;
    case Dir::W: *bdy = Dir::NW; break;
    default: assert(false);
    }
    *y -= d.dy;
  } else if (*y < d.y) {
    if (_b.GetBC(Dir::S) != Boundaries::bc_periodic) return false;
    switch (*bdy) {
    case Dir::inside: *bdy = Dir::S; break;
    case Dir::E: *bdy = Dir::SE; break;
    case Dir::W: *bdy = Dir::SW; break;
    default: assert(false);
    }
    *y += d.dy;
  }
  assert(InsideRect(d, *x, *y));
  return true;
}

// Move in direction bdy.
inline void RmuAnalyzer::
MapToPeriodic (Dir::Enum bdy, double* x, double* y) const {
  if (bdy == Dir::inside) return;
  const Rect& d = _rmu.GetDomain();
  switch (bdy) {
  case Dir::SE: case Dir::E: case Dir::NE: *x += d.dx; break;
  case Dir::NW: case Dir::W: case Dir::SW: *x -= d.dx; break;
  default: break;
  }
  switch (bdy) {
  case Dir::NE: case Dir::N: case Dir::NW: *y += d.dy; break;
  case Dir::SW: case Dir::S: case Dir::SE: *y -= d.dy; break;
  default: break;
  }
}

// Move in direction bdy.
inline void RmuAnalyzer::MapToPeriodic (Dir::Enum bdy, int* x, int* y) const {
  if (bdy == Dir::inside) return;
  const Rect& d = _rmu.GetDomain();
  switch (bdy) {
  case Dir::SE: case Dir::E: case Dir::NE: *x += UnitDx(d.dx); break;
  case Dir::NW: case Dir::W: case Dir::SW: *x -= UnitDx(d.dx); break;
  default: break;
  }
  switch (bdy) {
  case Dir::NE: case Dir::N: case Dir::NW: *y += UnitDy(d.dy); break;
  case Dir::SW: case Dir::S: case Dir::SE: *y -= UnitDy(d.dy); break;
  default: break;
  }
}

void RmuAnalyzer::GetPeriodicTri2 (const TriTag& tag, double* v1, double* v2,
                                   double* v3) const {
  assert(tag.id > 0 && tag.id <= _tri.Size(2));
  double* vs[3];
  vs[0] = v1; vs[1] = v2; vs[2] = v3;
  for (size_t i = 0; i < 3; i++) {
    if (!vs[i]) continue;
    MapToPeriodicWithTriAnchoredAt(tag.anchor, i, &_tri(1, tag.id),
                                   &_tri_bdy(1, tag.id), vs[i], vs[i]+1);
    MapToPeriodic(tag.dir, vs[i], vs[i]+1);
  }
}

const Matrix<Dir::Enum>& RmuAnalyzer::GetTriBdy () const { return _tri_bdy; }

// Map tri so that tri[anchor_id] is Dir::inside the domain. Return the mapped
// (x,y) of vertex v_id.
inline void RmuAnalyzer::MapToPeriodicWithTriAnchoredAt (
  size_t anchor_id, size_t v_id, const RectId tri[3],
  const Dir::Enum tri_bdy[3], double* x, double* y) const
{
  *x = _x(1, tri[v_id]+1);
  *y = _x(2, tri[v_id]+1);
  MapToPeriodic(tri_bdy[v_id], x, y);
  MapToPeriodic(Dir::Opposite(tri_bdy[anchor_id]), x, y);
}

Dir::Enum DirectionToXY (const Rect& r, double x, double y) {
  Dir::Enum dir1 = Dir::inside, dir2 = Dir::inside;
  if (x < r.x) dir1 = Dir::W;
  else if (x > r.x + r.dx) dir1 = Dir::E;
  if (y < r.y) dir2 = Dir::S;
  else if (y > r.y + r.dy) dir2 = Dir::N;
  return Dir::Combine(dir1, dir2);
}

// Set *r to the periodic rectangle that covers (x, y). Return false if
// impossible.
bool RmuAnalyzer::MapToCover (double x, double y, Rect* r) const {
  Dir::Enum dir = DirectionToXY(*r, x, y);
  if (dir == Dir::inside) return true;
  Dir::Enum dir1, dir2;
  Dir::BreakUp(dir, dir1, dir2);
  if (dir1 != Dir::inside && _b.GetBC(dir1) != Boundaries::bc_periodic)
    return false;
  if (dir2 != Dir::inside && _b.GetBC(dir2) != Boundaries::bc_periodic)
    return false;
  switch (dir1) {
  case Dir::S: r->y -= _r.dy; break;
  case Dir::N: r->y += _r.dy; break;
  default: break;
  }
  switch (dir2) {
  case Dir::E: r->x += _r.dx; break;
  case Dir::W: r->x -= _r.dx; break;
  default: break;
  }
  return InsideRect(*r, x, y);
}

void RmuAnalyzer::Distance (RectId id1, RectId id2, double* dx, double* dy,
                            Dir::Enum* dir) const {
  assert(id1 < _nr || id2 < _nr);
  if (dir) *dir = Dir::inside;
  int Dx = 0, Dy = 0;
  if (id1 < _nr && id2 < _nr) {
    const Rect* r1 = &_rmu.GetRects()[id1];
    const Rect* r2 = &_rmu.GetRects()[id2];
    PointRect::Distance(PointRect(*r1), PointRect(*r2), dx, dy);
    if (NoPeriodic(_b)) return;
    if (_b.GetBC(Dir::E) == Boundaries::bc_periodic) {
      int sign = 1;
      if (r2->x < r1->x) {
        std::swap(r1, r2);
        sign = -1;
      }
      double tdx;
      PointRect::Distance(PointRect(*r1),
                          PointRect(*r2, -_rmu.GetDomain().dx, 0),
                          &tdx, NULL);
      if (tdx < *dx) {
        *dx = tdx;
        Dx = -sign;
      }
      if (sign < 0) std::swap(r1, r2);
    }
    if (_b.GetBC(Dir::N) == Boundaries::bc_periodic) {
      int sign = 1;
      if (r2->y < r1->y) {
        std::swap(r1, r2);
        sign = -1;
      }
      // Don't have to include delta_x since we get distances separately in each
      // dimension.
      double tdy;
      PointRect::Distance(PointRect(*r1),
                          PointRect(*r2, 0, -_rmu.GetDomain().dy),
                          NULL, &tdy);
      if (tdy < *dy) {
        *dy = tdy;
        Dy = -sign;
      }
    }
  } else {
    int sign = 1;
    if (id2 < _nr) {
      std::swap(id1, id2);
      sign = -1;
    }
    const Rect& r = _rmu.GetRects()[id1];
    const PointRect pr(r);
    Point p(_x(1, id2+1), _x(2, id2+1));
    PointRect::Distance(pr, p, dx, dy);
    if (NoPeriodic(_b)) return;
    // We're here only if there are edge points, which in turn implies there
    // is at least one nonperiodic boundary.
    assert(_b.GetBC(Dir::E) != Boundaries::bc_periodic ||
           _b.GetBC(Dir::N) != Boundaries::bc_periodic);
    if (_b.GetBC(Dir::E) == Boundaries::bc_periodic) {
      int signx = r.x < p.x ? -1 : 1;
      const double delta_x = signx * _rmu.GetDomain().dx;
      double tdx;
      PointRect::Distance(r, Point(p.x + delta_x, p.y), &tdx, NULL);
      if (tdx < *dx) {
        *dx = tdx;
        Dx = sign * signx;
      }
    } else {
      assert(_b.GetBC(Dir::N) == Boundaries::bc_periodic);
      int signy = r.y < p.y ? -1 : 1;
      const double delta_y = signy * _rmu.GetDomain().dy;
      double tdy;
      PointRect::Distance(r, Point(p.x, p.y + delta_y), NULL, &tdy);
      if (tdy < *dy) {
        *dy = tdy;
        Dy = sign * signy;
      }
    }
  }
  if (dir) *dir = Dir::ToEnum(Dx, Dy);
}

namespace {
  inline Rect& TranslateRect (const Rect& domain, const RectTag& rt, Rect& r) {
    r.x += domain.dx * rt.Dx;
    r.y += domain.dy * rt.Dy;
    return r;
  }

  inline Point& TranslatePoint (const Rect& domain, const RectTag& rt, Point& p)
  {
    p.x += domain.dx * rt.Dx;
    p.y += domain.dy * rt.Dy;
    return p;    
  }
}

void RmuAnalyzer::
Distance (const RectTag& rt1, const RectTag& rt2, double* dx, double* dy) const
{
  assert(rt1.id < _nr || rt2.id < _nr);
  if (rt1.id < _nr && rt2.id < _nr) {
    Rect r1 = _rmu.GetRects()[rt1.id];
    TranslateRect(_rmu.GetDomain(), rt1, r1);
    Rect r2 = _rmu.GetRects()[rt2.id];
    TranslateRect(_rmu.GetDomain(), rt2, r2);
    PointRect::Distance(PointRect(r1), PointRect(r2), dx, dy);
  } else {
    if (rt2.id < _nr) return Distance(rt2, rt1, dx, dy);
    Rect r1 = _rmu.GetRects()[rt1.id];
    TranslateRect(_rmu.GetDomain(), rt1, r1);
    Point p(_x(1, rt2.id+1), _x(2, rt2.id+1));
    TranslatePoint(_rmu.GetDomain(), rt2, p);
    PointRect::Distance(PointRect(r1), p, dx, dy);
  }
}

namespace bdy {
struct BdyPoint {
  int x, y;
  RectId id;
  BdyPoint() : x(0), y(0), id(-1) {}
  BdyPoint (int ix, int iy, RectId iid) { Set(ix, iy, iid); }
  void Set (int ix, int iy, RectId iid) { x = ix; y = iy; id = iid; }
};

bool operator< (const BdyPoint& ep1, const BdyPoint& ep2) {
  if (ep1.x < ep2.x) return true;
  else if (ep1.x > ep2.x) return false;
  else return (ep1.y < ep2.y);
}

bool AddBdyPointIfNew (const RmuAnalyzer& ra, set<BdyPoint>* eps, double x,
                       double y, RectId* avlbl_id, RectId* id) {
  Dir::Enum bdy = ra.OnWhichBoundary(x, y);
  assert(bdy != Dir::inside);
  if (!ra.GetBoundaries().HasVBC(bdy)) return false;

  BdyPoint ep(ra.UnitX(x), ra.UnitY(y), -1);
  set<BdyPoint>::const_iterator it = eps->find(ep);
  if (it == eps->end()) {
    ep.id = *avlbl_id;
    eps->insert(ep);
    (*avlbl_id)++;
    *id = ep.id;
  } else {
    *id = it->id;
  }
  return true;
}

// Gather ye boundary points that are on v(0)-BC boundaries.
inline size_t
GatherBdyPoints (const RmuAnalyzer& ra, set<BdyPoint>* eps, Dir::Enum dir,
                 const Rect&r, RectId* next_id, RectId ids[3],
                 Dir::Enum bdys[3]) {
#define aepin(x, y, bdy)                                        \
  if (AddBdyPointIfNew(ra, eps, x, y, next_id, ids + nid)) {    \
    bdys[nid] = bdy;                                            \
    nid++;                                                      \
  }
  const double rdx = r.x + r.dx;
  const double rdy = r.y + r.dy;
  size_t nid = 0;
  switch (dir) {
  case Dir::E:  aepin(rdx, rdy, Dir::E); aepin(rdx, r.y, Dir::E); break;
  case Dir::NE: aepin(rdx, r.y, Dir::E); aepin(rdx, rdy, Dir::NE);
    aepin(r.x, rdy, Dir::N); break;
  case Dir::N:  aepin(rdx, rdy, Dir::N); aepin(r.x, rdy, Dir::N); break;
  case Dir::NW: aepin(rdx, rdy, Dir::N); aepin(r.x, rdy, Dir::NW);
    aepin(r.x, r.y, Dir::W); break;
  case Dir::W:  aepin(r.x, rdy, Dir::W); aepin(r.x, r.y, Dir::W); break;
  case Dir::SW: aepin(r.x, rdy, Dir::W); aepin(r.x, r.y, Dir::SW);
    aepin(rdx, r.y, Dir::S); break;
  case Dir::S:  aepin(r.x, r.y, Dir::S); aepin(rdx, r.y, Dir::S); break;
  case Dir::SE: aepin(rdx, rdy, Dir::E); aepin(rdx, r.y, Dir::SE);
    aepin(r.x, r.y, Dir::S); break;
  case Dir::inside: break;
  }
  return nid;

#undef aepin
#undef rdx
#undef rdy
}
} // bdy

// Build the adjacency graph over rectangles and boundary points.
void RmuAnalyzer::OrganizeInSpace (Vve& age) {
  using namespace bdy;
  assert(_nr > 0);

  age.resize(_nr);
  _x.Resize(2, _nr);
  set<BdyPoint> bdy_pts;
  RectId next_id = _nr;
  vector<bool> handled_corner(4, false);
  for (RectId id = 0; id < _nr; id++) {
    // Get this rectangle's neighbors.
    GetNbrs(id, age[id]);

    // Connect rectangle centers.
    const Rect& r = _rs[id];
    _x(1, id+1) = r.x + 0.5*r.dx;
    _x(2, id+1) = r.y + 0.5*r.dy;

    // Handle the boundary conditions. If a rectangle is on a bc_(0)velocity
    // boundary, add the boundary-touching corner points of the rectangle to
    // the graph.
    //   First, see if we can can continue.
    Dir::Enum bdy = OnWhichBoundary(r);
    if (bdy == Dir::inside || !_b.HasVBC(bdy)) continue;
    //   We're still here, so that means both of the 2 (cardinal edge) or at
    // least 1 of the 3 (corner) points is on a v-BC boundary.
    RectId ids[3];
    Dir::Enum bdys[3];
    //    Next, gather bdy points that are on a v-BC boundary.
    size_t nid = GatherBdyPoints(*this, &bdy_pts, bdy, r, &next_id, ids, bdys);
    assert(nid >= 2);
    //   Finally, add the rectangle-boundary point and bdy point-bdy point
    // edges. This is complicated because we need to handle several cases:
    // rectangle on interior v(0)-BC boundary; rectangle at corner of
    // all-v(0)-BC boundaries; rect at corner with one v(0)-BC and one
    // periodic. For the two corner cases, we need to distinguish between the
    // boundary point that is itself on the corner and the 2 (all v(0)-BC) or
    // 1 (mix) non-corner points.
    for (size_t i = 0; i < nid; i++) {
      RectId idi = ids[i];
      // Connect rect center and this edge point.
      age[id].push_back(Edge(idi));
      assert(idi <= age.size());
      if (idi == age.size()) {
        age.push_back(vector<Edge>());
        _bdypt_side.push_back(Dir::inside);
      }
      _bdypt_side[idi - _nr] = bdys[i];
      age[idi].push_back(Edge(id));
      // Connect edge points to each other.
      if (Dir::IsCardinal(bdy)) {
        assert(nid == 2);
        age[idi].push_back(Edge(ids[(i + 1) % 2]));
      } else { // Rect r is on a corner.
        if (Dir::IsCorner(bdys[i])) { // This edge point is on the corner.
          for (size_t j = 0; j < nid; j++)
            if (j != i) age[idi].push_back(Edge(ids[j]));
          if (!NoPeriodic(_b)) // One side is periodic and the other is not.
            OrganizeInSpaceHandleBcCornerVP(id, idi, bdy, handled_corner, age);
        } else { // This edge point is not on the corner.
          for (size_t j = 0; j < nid; j++)
            if (Dir::IsCorner(bdys[j])) {
              age[idi].push_back(Edge(ids[j]));
              break;
            }
        }
      }
    } // corner point i loop
  } // rect id loop

  // Add the boundary points to _x.
  _x.Reshape(2, _nr + bdy_pts.size());
  for (set<BdyPoint>::const_iterator it = bdy_pts.begin();
       it != bdy_pts.end(); ++it) {
    _x(1, it->id+1) = DomainX(it->x);
    _x(2, it->id+1) = DomainY(it->y);
  }

  // Sort for fast finding later on.
  for (RectId id = 0; id < age.size(); id++)
    std::sort(age[id].begin(), age[id].end());
}

// Handle a corner rectangle with one side (0)v-BC and the other periodic.
// Corner points when two sides are periodic and two are (0)v-BC are connected
// to periodic rectangle centers.
void RmuAnalyzer::OrganizeInSpaceHandleBcCornerVP (
  const RectId rect_id, // the current rectangle
  const RectId node_id, // the vertex of the rectangle
  const Dir::Enum bdy,  // the corner the rectangle is on
  vector<bool>& handled_corner, // avoid duplication that GetEdges can't detect
  Vve& age)
{
  Dir::Enum bdy1, bdy2;
  Dir::BreakUp(bdy, bdy1, bdy2);
  Dir::Enum pbdy, vbdy;
  if (_b.GetBC(bdy1) == Boundaries::bc_periodic) {
    pbdy = bdy1; vbdy = bdy2;
  } else {
    pbdy = bdy2; vbdy = bdy1;
  }
  // Find the rectangle at the opposite corner.
  int pid = -1;
  for (size_t k = 0; k < age[rect_id].size(); k++) {
    RectId kid = age[rect_id][k].id;
    if (kid >= _nr) continue;
    if (age[rect_id][k].bdy != Dir::inside &&
        !Dir::IsCardinal(OnWhichBoundary(_rs[kid]))) {
      pid = kid;
      break;
    }}
  assert(pid >= 0 && Dir::IsCardinal(pbdy));
  // Make the edges.
  age[pid].push_back(Edge(node_id, Dir::Opposite(pbdy)));
  age[node_id].push_back(Edge(pid, pbdy));
  if (!handled_corner[(int) vbdy])
    handled_corner[(int) vbdy] = true;
  else {
    age[pid].back().use_in_GetEdge = false;
    age[node_id].back().use_in_GetEdge = false;
  }
}

class PerimeterCrawler {
public:
  // (ux, uy) is the unit cell size.
  PerimeterCrawler(const Rect* r, double ux, double uy);

  // Get the crawler's current point.
  double cx () const { return _pci.cx(); }
  double cy () const { return _pci.cy(); }
  // Returns false without moving the crawler when it's made it once around
  // qt's rect.
  bool Incr();
  // Cover part of the perimeter with this rect.
  void Cover(const Rect& r);

private:
  class PcIncr {
  public:
    PcIncr(const Rect* r, double ux, double uy);

    const Rect& GetRect () const { return _r; }
    // Returns false if either the crawler is back to the beginning or the
    // whole perimeter is covered.
    bool Incr();
    double cx () const { return _cx; }
    double cy () const { return _cy; }
    double ux () const { return _ux; }
    double uy () const { return _uy; }

  private:
    enum Region { E = 0, NE, N, NW, W, SW, S, SE };
    const Rect& _r;
    double _ux, _uy;
    double _cx, _cy;
    Region _region;
  };

private:
  PcIncr _pci;
  bool _done;
  vector<bool> _is_covered;
  size_t _n_covered, _cover_i;
};

PerimeterCrawler::PcIncr::PcIncr (const Rect* r, double ux, double uy)
  : _r(*r), _ux(ux), _uy(uy)
{
  _cx = _r.x + _r.dx + 0.5*_ux;
  _cy = _r.y - 0.5*_uy;
  _region = SE;
}

bool PerimeterCrawler::PcIncr::Incr () {
  // A unit cell is 2x2.
  switch (_region) {
  case SE: _cy += _uy; _region = E; break;
  case E:  _cy += _uy; if (_cy > _r.y + _r.dy) _region = NE; break;
  case NE: _cx -= _ux; _region = N; break;
  case N:  _cx -= _ux; if (_cx < _r.x) _region = NW; break;
  case NW: _cy -= _uy; _region = W; break;
  case W:  _cy -= _uy; if (_cy < _r.y) _region = SW; break;
  case SW: _cx += _ux; _region = S; break;
  case S:  _cx += _ux; if (_cx > _r.x + _r.dx) _region = SE; break;
  }
  return _region != SE;
}

PerimeterCrawler::PerimeterCrawler (const Rect* r, double ux, double uy)
  : _pci(r, ux, uy)
{
  _is_covered.resize((size_t) (2*(r->dx/ux + r->dy/uy) + 4), false);
  _n_covered = 0;
  _cover_i = 0;
}

bool PerimeterCrawler::Incr () {
  if (_n_covered == _is_covered.size() ||
      _cover_i == _is_covered.size()) return false;
  do {
    _cover_i++;
    if (_cover_i == _is_covered.size()) return false;
    _pci.Incr();
  } while (_is_covered[_cover_i]);
  return true;
}

//assume intersect(_pci.GetRect(), r) is empty.
void PerimeterCrawler::Cover (const Rect& r) {
  //todo-opt
  PcIncr pci(&_pci.GetRect(), _pci.ux(), _pci.uy());
  for (size_t k = 0; ; k++) {
    if (InsideRect(r, pci.cx(), pci.cy())) {
      if (!_is_covered[k]) _n_covered++;
      _is_covered[k] = true;
    }
    if (!pci.Incr()) break;
  }
}

void RmuAnalyzer::GetNbrs (RectId id, vector<Edge>& nbrs) const {
  //todo-opt Take a non-empty nbrs and fill pc with what we already know.
  PerimeterCrawler pc(&_rs[id], _ux, _uy);
  do {
    //todo-opt Replace with a method that goes up the tree and then down. This
    // one always starts at the top.
    int id = _rmu.GetRectId(pc.cx(), pc.cy());
    if (id >= 0) {
      // pc's current point is in the domain.
      nbrs.push_back(Edge(id));
      pc.Cover(_rs[id]);
    } else {
      // pc's current point is outside the domain. We need to handle it if we
      // have periodic BC here.
      Dir::Enum bdy;
      double x = pc.cx(), y = pc.cy();
      if (MapToDomain(&x, &y, &bdy)) {
        // Moving in the direction bdy brings us to the new point (x, y).
        id = _rmu.GetRectId(x, y);
        assert(id >= 0);
        nbrs.push_back(Edge(id, bdy));
        Rect cr = _rs[id];
#ifndef NDEBUG
        bool ret =
#endif
          MapToCover(pc.cx(), pc.cy(), &cr);
        assert(ret);
        pc.Cover(cr);
      }
    }
  } while (pc.Incr());
}

const Boundaries& RmuAnalyzer::GetBoundaries () const { return _b; }

const vector<RectId>& RmuAnalyzer::GetNbrs (RectId id) const {
  assert(id < _ag.size());
  return _ag[id];
}
const vector<Dir::Enum>& RmuAnalyzer::GetNbrsBdys (RectId id) const {
  assert(id < _ag.size());
  return _ag_bdy[id];
}

Dir::Enum RmuAnalyzer::GetBoundaryDir (RectId id) const {
  if (id < _nr) return Dir::inside;
  assert(_bdypt_side[id - _nr] != Dir::inside);
  return _bdypt_side[id - _nr];
}

namespace {
void PrintEdges (RectId id, Dir::Enum corner, const vector<Edge>& es) {
  printf("  point %d at corner %d: [", (int) id, (int) corner);
  for (size_t i = 0; i < es.size(); i++)
    printf("(%d %d %d) ", (int) es[i].id, (int) es[i].bdy,
           (int) es[i].use_in_GetEdge);
  printf("];\n");
}
}

// -----------------------------------
// Triangulation.

class Triangulator {
public:
  Triangulator (const RmuAnalyzer& irmu, const Vve& rage,
                const Matrix<double>& x, const Boundaries& b,
                double domain_size)
    : x(x), rmu(irmu), _rage(rage), _b(b),
      _atol(1e2 * domain_size * std::numeric_limits<double>::epsilon())
  {}

  void Triangulate(Matrix<Edge>* trie);
  double AspectRatio(const Edge trie[3], size_t* max_edge_len_idx) const;
  void MakeCcw(Edge a[3]) const;

  struct RectData {
    RectId id;
    double sz; // RectU size

    void Set (RectId iid, double dx) { id = iid; sz = dx; }
    bool operator< (const RectData& r2) const { return sz < r2.sz; }
  };
  vector<RectData> rd;
  const Matrix<double>& x;  // [rect centers, edge points]
  const RmuAnalyzer& rmu;

private:
  const Vve& _rage; // RectId adjacency graph.
  const Boundaries& _b;
  double _atol;

private:
  void GetEdges(Vve& age);
  void MakeTri(const Vve& age, Matrix<Edge>& trie);
  bool IntersectWithPrevious(const Vve& age, RectId iid, const Edge& je)
    const;
};

void Triangulator::Triangulate (Matrix<Edge>* trie) {
  Vve age;
  Timer* t = Ca::GetTimer();
  t->Tic(0);
  GetEdges(age);
  double et_ge = t->Toc(0); t->Tic(0);
  MakeTri(age, *trie);
  double et_mt = t->Toc(0);
  caprint("TRI: ge (%1.2e) 1.00 mt %1.2f\n", et_ge, et_mt/et_ge);
}

// From the adjacency graph, extract nonintersecting edges to form triangles.
void Triangulator::GetEdges (Vve& age) {
  size_t nr = this->rd.size();
  size_t nx = this->x.Size(2);
  assert(nx >= nr);

  // Set up rd and iperm.
  vector<size_t> iperm(nx);
  // Sort rects in increasing size.
  std::sort(this->rd.begin(), this->rd.end());
  // Inverse permutation.
  for (size_t i = 0; i < nr; i++) iperm[this->rd[i].id] = i;
  // Add edge points.
  this->rd.insert(this->rd.end(), nx - nr, Triangulator::RectData());
  for (size_t i = nr; i < nx; i++) {
    this->rd[i].id = (RectId) i;
    // Don't need to set: this->rd[i].sz
    iperm[i] = i;
  }

  // Test if line segment i-j intersects line segs already established.
  age.resize(nx);
  for (size_t i = 0; i < nx; i++) {
    RectId iid = this->rd[i].id;
    const vector<Edge>& inbrs = _rage[iid];
    for (size_t j = 0; j < inbrs.size(); j++) {
      if (!inbrs[j].use_in_GetEdge) continue;
      const Edge je(inbrs[j].id, inbrs[j].bdy);
      const size_t jip = iperm[je.id];
      // Connect to only those nodes that have been processed. This has two
      // effects: 1. It prevents both directional edges between two nodes from
      // being used, leading to redundant triangles. (We could write
      // IntersectWithPrevious to account for these, but why bother when an
      // ordering gives it to us for free.) 2. It goes from small to large,
      // which is a nice way of getting mostly high-quality triangles.
      // (FlipEdges, if we use it, takes care of the rest.)
      if (i < jip) continue;
      assert(i != jip);
      if (!IntersectWithPrevious(age, iid, je)) age[iid].push_back(je);
    }
  }
}

//todo-numerical Eventually, I want to transition all the calculations to unit
// coords to avoid roundoff. For now, I convert here.
bool Triangulator::
IntersectWithPrevious (const Vve& age, RectId iid, const Edge& je) const {
  // Proposed line segment. Point i is inside the domain.
  PointT<int> jp(rmu.UnitX(this->x(1, je.id+1)),
                 rmu.UnitY(this->x(2, je.id+1)));
  rmu.MapToPeriodic(je.bdy, &jp.x, &jp.y);
  PointLineSegT<int> lsij(PointT<int>(rmu.UnitX(this->x(1, iid+1)),
                                      rmu.UnitY(this->x(2, iid+1))),
                          jp);
  for (size_t in = 0; in < _rage[iid].size(); in++) {
    const RectId inid = _rage[iid][in].id;
    PointT<int> inp(rmu.UnitX(this->x(1, inid+1)),
                    rmu.UnitY(this->x(2, inid+1)));
    rmu.MapToPeriodic(_rage[iid][in].bdy, &inp.x, &inp.y);
    for (size_t k = 0; k < age[inid].size(); k++) {
      const RectId kid = age[inid][k].id;
      // Does the proposed line seg (iid-je.id) intersect an established one
      // (inid-kid)?
      PointT<int> kp(rmu.UnitX(this->x(1, kid+1)),
                     rmu.UnitY(this->x(2, kid+1)));
      rmu.MapToPeriodic(age[inid][k].bdy, &kp.x, &kp.y);
      rmu.MapToPeriodic(_rage[iid][in].bdy, &kp.x, &kp.y);
      PointLineSegT<int> lsink(inp, kp);
      if (PointLineSeg::IntersectOpen(lsij, lsink)) return true;
    }
  }
  return false;
}

namespace maketri {
  inline bool IsIn (const vector<Edge>& set, RectId id, Dir::Enum* bdy) {
    //todo-opt Sort ag[i], ag_bdy[i] for each i when constructed, and then do
    // binary search here.
    for (size_t i = 0; i < set.size(); i++)
      if (set[i].id == id) {
        *bdy = set[i].bdy;
        return true;
      }
    return false;
  }
}

// Turn the set of edges into individual triangles.
void Triangulator::MakeTri (const Vve& age, Matrix<Edge>& trie) {
  vector<Edge> vtrie;
  size_t k = 0;
  const double* x = this->x.GetPtr();
  for (RectId n1 = 0; n1 < age.size(); n1++) {
    for (int in2 = 0; in2 < (int) age[n1].size() - 1; in2++) {
      const RectId n2 = age[n1][in2].id;
      for (size_t in3 = in2 + 1; in3 < age[n1].size(); in3++) {
        const RectId n3 = age[n1][in3].id;
        Dir::Enum bdy;
        if (maketri::IsIn(age[n2], n3, &bdy) ||
            maketri::IsIn(age[n3], n2, &bdy)) {
          // Convention: We pin vertex 1 to be inside the domain.
          vtrie.push_back(Edge(n1, Dir::inside));
          vtrie.push_back(Edge(n2, age[n1][in2].bdy));
          vtrie.push_back(Edge(n3, age[n1][in3].bdy));
          k += 3;
          // Triangles are built using _x. Let's call the original triangle
          // defined simply on x as the "base tri". Then _tri_bdy tells us
          // whether certain corners should be mapped periodically. If only
          // one boundary is involved, the base tri is flipped once, thus
          // reversing its vertex ordering. If two boundaries (the most
          // possible) are involved, then it's flipped twice, returning it to
          // its original vertex ordering. Here I handle these details.
          const bool all_inside = age[n1][in2].bdy == Dir::inside
            && age[n1][in3].bdy == Dir::inside;
          bool uses_2bdys = !all_inside
            && age[n1][in2].bdy != Dir::inside
            && age[n1][in3].bdy != Dir::inside
            && age[n1][in2].bdy != age[n1][in3].bdy;
          assert(!all_inside || bdy == Dir::inside);
          const bool is_ccw = Tri2::IsCcw(x + 2*n1, x + 2*n2, x + 2*n3);
          if (( (all_inside || uses_2bdys) && !is_ccw) ||
              (!(all_inside || uses_2bdys) &&  is_ccw))
            std::swap(vtrie[k - 2], vtrie[k - 1]);
        }}}
  }
  trie.Resize(3, k/3);
  memcpy(trie.GetPtr(), &vtrie[0], k*sizeof(Edge));
}

inline void Triangulator::MakeCcw (Edge a[3]) const {
  Tri2 t2(&x(1, a[0].id+1), &x(1, a[1].id+1), &x(1, a[2].id+1));
  const bool is_ccw = t2.IsCcw();
  const bool all_inside = a[0].bdy == Dir::inside && a[1].bdy == Dir::inside
    && a[2].bdy == Dir::inside;
  bool uses_2bdys = false;
  if (!all_inside) {
    // See if the triangle goes through two boundaries.
    Dir::Enum bdy1 = Dir::inside;
    for (size_t i = 0; i < 3; i++) {
      if (a[i].bdy != Dir::inside) {
        if (bdy1 == Dir::inside)
          bdy1 = a[i].bdy;
        else {
          uses_2bdys = a[i].bdy != bdy1;
          break;
        }}}
  }
  if (( (all_inside || uses_2bdys) && !is_ccw) ||
      (!(all_inside || uses_2bdys) &&  is_ccw))
    std::swap(a[1], a[2]);
}

double Triangulator::
AspectRatio (const Edge trie[3], size_t* max_edge_len_idx) const {
  double vs[6];
  for (size_t i = 0; i < 3; i++) {
    memcpy(vs + 2*i, &x(1, trie[i].id+1), 2*sizeof(double));
    rmu.MapToPeriodic(trie[i].bdy, vs + 2*i, vs + 2*i+1);
  }
  Tri2 t2(vs, vs+2, vs+4);
  return t2.AspectRatio(max_edge_len_idx);
}

void RmuAnalyzer::Triangulate (const Vve& age, Matrix<Edge>& trie) {
  assert(_nr > 0 && age.size() > 0 && _x.Size(2) >= _nr);
  Triangulator t(*this, age, _x, _b, std::max(_r.dx, _r.dy));
  t.rd.resize(_nr);
  for (size_t i = 0; i < _nr; i++) {
    const Rect& r = _rs[i];
    t.rd[i].Set(i, r.dx);
  }
  t.Triangulate(&trie);
}

namespace r2t {
typedef vector< vector<size_t> > Vvi;

// On output, trie(:, tis[i]) are the triangles associated with _x(:,i) for i
// <= _nr.
void DistributeTris (const Matrix<RectId>& tri,
                     const Matrix<Dir::Enum>& tri_bdy, Vvi& tis) {
  size_t nr = tis.size();
  assert(nr > 0);
  const RectId* ptri = tri.GetPtr();
  for (size_t i = 0, n = tri.Size(); i < n; i++)
    if (ptri[i] < nr) tis[ptri[i]].push_back(i / 3);
}
}

// Angle of the trailing arm in a CCW direction, anchored at rect center ir.
double RmuAnalyzer::
Angle (RectId ir, const RectId tri[3], const Dir::Enum tri_bdy[3]) const {
  // Get the triangle vertex at rect ir's center.
  size_t iv;
  for (iv = 0; iv < 3; iv++) if (tri[iv] == ir) break;
  assert(tri[iv] == ir);
  // Trailing edge goes from iv to jv.
  const size_t jv = (iv + 1) % 3;

  double x2[2];
  memcpy(x2, &_x(1, tri[jv]+1), 2*sizeof(double));
  MapToPeriodicWithTriAnchoredAt(iv, jv, tri, tri_bdy, x2, x2+1);

  const double dx = x2[0] - _x(1, tri[iv]+1);
  const double dy = x2[1] - _x(2, tri[iv]+1);
  return atan2(dy, dx);
}

bool RmuAnalyzer::RtoT::operator< (const RtoT& r1) const
{ return theta < r1.theta; }

// Fill in _r2t. This data structure allows fast query of the triangle that
// covers a point. First the covering rectangle is found in O(log n)
// time. Then the triangle is found in O(1) steps using _r2t.
void RmuAnalyzer::MakeRtoT () {
  assert(_nr > 0 && _tri.Size() > 0 && _tri.Size() == _tri_bdy.Size());

  // For each rect, get a list of triangles with a vertex at the rect center.
  r2t::Vvi tis(_nr);
  r2t::DistributeTris(_tri, _tri_bdy, tis);

  _r2t.resize(_nr);
  for (size_t ir = 0; ir < _nr; ir++) { // For each rect
    size_t nt = tis[ir].size();
    assert(nt > 0);
    vector<RtoT>& r2ts = _r2t[ir];
    r2ts.resize(nt);
    for (size_t it = 0; it < nt; it++) { // For each triangle in the rect
      size_t ti = tis[ir][it];
      r2ts[it].ti = ti; // Triangle index
      // Tri arm angle.
      r2ts[it].theta = Angle(ir, &_tri(1, ti+1), &_tri_bdy(1, ti+1));
    }
    std::sort(r2ts.begin(), r2ts.end()); // Ascending arm angle

#ifndef NDEBUG
    // MakeRtoT does complicated stuff with all the tris and rects, so it's the
    // best and last place to catch problems. Let's keep this test and thorough
    // debug output in the code indefinitely.
    for (size_t it = 0; it < nt-1; it++) {
      if (r2ts[it].theta == r2ts[it+1].theta) {
        // Ruh roh!
        printf("ir = %d (%f, %f)\n", (int) ir, _x(1, ir+1), _x(2, ir+1));
        printf("tri = [");
        for (size_t it = 0; it < nt; it++) {
          size_t ti = r2ts[it].ti;
          printf("%d %d %d\n",
                 (int) _tri(1, ti+1), (int) _tri(2, ti+1), (int) _tri(3, ti+1));
        }
        printf("];\n");
        for (size_t it = 0; it < nt; it++) {
          size_t iv;
          size_t ti = r2ts[it].ti;
          for (iv = 0; iv < 3; iv++) if (_tri(iv+1, ti+1) == ir) break;
          assert(_tri(iv+1, ti+1) == ir);
          Tri2 t2(&_x(1, _tri(1, ti+1)+1), &_x(1, _tri(2, ti+1)+1),
                  &_x(1, _tri(3, ti+1)+1));
          printf("%5.1f | %d %d %d | ccw %d \n", r2ts[it].theta/M_PI*180,
                 (int) _tri_bdy(iv+1, ti+1),
                 (int) _tri_bdy(((iv+1)%3)+1, ti+1),
                 (int) _tri_bdy(((iv+2)%3)+1, ti+1),
                 t2.IsCcw());
        }
        printf("\n");
        assert(false);
      }
    }
#endif
#ifndef NDEBUG
    for (size_t it = 0; it < nt-1; it++) assert(r2ts[it] < r2ts[it+1]);
#endif
  }
}

// If it exists, return the triangle that contains (x, y). tag contains the data
// necessary to translate the base tri.
bool RmuAnalyzer::GetTriTag (double x, double y, TriTag* tag) const {
  assert(tag);
  //todo-periodic-complete If I ever decide to allow neighborhoods to grow
  // beyond the size of one domain, then I need to generalize from dir to a
  // vector to allow multiple domain translations.
#ifndef NDEBUG
  bool ret =
#endif
    MapToDomain(&x, &y, &tag->dir);
  assert(ret);
  tag->dir = tag->dir;
  // Get the rect.
  const int rid = _rmu.GetRectId(x, y);
  if (rid == -1) return false;
  assert(rid >= 0 && rid < (int) _r2t.size());
  // Get the triangle.
  const double dx = x - _x(1, rid+1);
  const double dy = y - _x(2, rid+1);
  RtoT key;
  key.theta = atan2(dy, dx);
  vector<RtoT>::const_iterator it =
    std::upper_bound(_r2t[rid].begin(), _r2t[rid].end(), key);
  if (it == _r2t[rid].end() || it == _r2t[rid].begin())
    tag->id = _r2t[rid].back().ti + 1;
  else {
    --it;
    tag->id = it->ti + 1;
  }
  // Get the triangle's anchor vertex.
  bool fnd = false;
  for (size_t i = 0; i < 3; i++)
    if ((int) _tri(i+1, tag->id) == rid) {
      tag->anchor = i;
      fnd = true;
      break;
    }
  assert(fnd);
  return true;
}

const Matrix<double>& RmuAnalyzer::GetX () const { return _x; }
const Matrix<RectId>& RmuAnalyzer::GetTri () const { return _tri; }

void RmuAnalyzer::GetTrisSharingRect (RectId id, vector<size_t>& tis) const {
  tis.clear();
  assert(id >= 0 && id < _r2t.size());
  for (size_t i = 0; i < _r2t[id].size(); i++)
    tis.push_back(_r2t[id][i].ti);
}

void RmuAnalyzer::
GetCornerBC (Dir::Enum corner, Dir::Enum& bdy, Boundaries::BC& bc) const {
  assert(Dir::IsCorner(corner));
  Dir::Enum dir1, dir2;
  Dir::BreakUp(corner, dir1, dir2);
  const Boundaries::BC bc1 = _b.GetBC(dir1);
  const Boundaries::BC bc2 = _b.GetBC(dir2);
#define do1 do { bc = bc1; bdy = dir1; } while (0)
#define do2 do { bc = bc2; bdy = dir2; } while (0)
  if      (bc1 == Boundaries::bc_velocity)  do1;
  else if (bc2 == Boundaries::bc_velocity)  do2;
  else if (bc1 == Boundaries::bc_0velocity) do1;
  else if (bc2 == Boundaries::bc_0velocity) do2;
  else if (bc1 == Boundaries::bc_free)      do1;
  else if (bc2 == Boundaries::bc_free)      do2;
  else assert(false);
#undef do1
#undef do2
}

//todo This procedure is not consistent with BrickMeshBemBuilder.cpp. However,
// I'm using it for external interpolation and not for IGA calculations, so I'm
// leaving it as is for now.
void RmuAnalyzer::
PrepInterp (vector<double>& z, const double bdy_vals[4]) const {
  size_t nx = _x.Size(2);
  if (_nr == nx) return;
  z.resize(nx, 0);
  for (size_t i = _nr; i < nx; i++) {
    Dir::Enum pdir = GetBoundaryDir(i);
    Dir::Enum side;
    Boundaries::BC bc;
    if (Dir::IsCorner(pdir))
      GetCornerBC(pdir, side, bc);
    else {
      side = pdir;
      bc = _b.GetBC(side);
    }
    switch (bc) {
    case Boundaries::bc_velocity:
      z[i] = bdy_vals[(int) side];
      break;
    case Boundaries::bc_free: {
      const vector<RectId>& ids = _ag[i];
      assert((Dir::IsCorner(pdir) && ids.size() == 3) || ids.size() == 4);
      double tot_area = 0;
      z[i] = 0;
      size_t cnt = 0;
      for (size_t j = 0; j < ids.size(); j++)
        if (ids[j] < _nr) {
          const Rect& r = _rs[ids[j]];
          const double area = r.dx * r.dy;
          tot_area += area;
          z[i] += area * z[ids[j]];
          cnt++;
        }
      assert(cnt == 2);
      z[i] /= tot_area;
    } break;
    default: break;
    }
  }
}

void RmuAnalyzer::
LinearInterp (const vector<double>& z,
              const vector<double>& xi, const vector<double>& yi,
              vector<double>& zi) const { 
  assert(z.size() == _x.Size(2));
  assert(xi.size() == yi.size());
  zi.resize(xi.size());
  for (size_t i = 0; i < xi.size(); i++) {
    TriTag tag;
    if (!GetTriTag(xi[i], yi[i], &tag)) {
      zi[i] = 0;
      continue;
    }
    const RectId* tri = &_tri(1, tag.id);
    double v1[2], v2[2], v3[2];
    GetPeriodicTri2(tag, v1, v2, v3);
    Tri2 t2(v1, v2, v3);
    double Ti[4];
    t2.BarycentricMatrix(Ti);
    double xy[2]; xy[0] = xi[i]; xy[1] = yi[i];
    double lambda[3];
    t2.ToBarycentric(Ti, xy, lambda);
    zi[i] = lambda[0]*z[tri[0]] + lambda[1]*z[tri[1]] + lambda[2]*z[tri[2]];
  }
}

void RmuAnalyzer::
InterpWExtrap (const vector<double>& z_int, const vector<double>& xi,
               const vector<double>& yi, vector<double>& zi) const {
  vector<double> z(z_int);
  size_t nx = _x.Size(2), nr = _rmu.GetRects().size();
  z.resize(nx);
  for (size_t i = nr; i < nx; i++) z[i] = Extrap(z_int, i);
  LinearInterp(z, xi, yi, zi);
}

void RmuAnalyzer::PrepExtrap (vector<double>& z) const {
  size_t nx = _x.Size(2), nr = _rmu.GetRects().size();
  z.resize(nx);
  for (size_t i = nr; i < nx; i++) z[i] = Extrap(z, i);
}

// Get the supporting points. This consists of two layers of elements starting
// with the one (corner) or two (side) adjacent to the boundary point. This
// procedure should give the same support as in BrickMeshBemBuilder for
// support=2.
static void
GetSupport (const vector< vector<RectId> >& ag, size_t nr, RectId id,
            vector<RectId>& sids) {
  set<RectId> support;
  for (size_t i = 0; i < ag[id].size(); i++) {
    RectId id1 = ag[id][i];
    if (id1 >= nr) continue;
    support.insert(id1);
    for (size_t j = 0; j < ag[id1].size(); j++)
      if (ag[id1][j] < nr) support.insert(ag[id1][j]);
  }

  const size_t n = support.size();
  assert(n >= 4);
  sids.reserve(n);
  for (set<RectId>::const_iterator it = support.begin(), end = support.end();
       it != end; ++it)
    sids.push_back(*it);
}

// Extrapolate the value at the bdy pt id.
double RmuAnalyzer::Extrap (const vector<double>& z_int, RectId id) const {
  vector<RectId> sids;
  GetSupport(_ag, _rmu.GetRects().size(), id, sids);
  const size_t n = sids.size();

  //todo-opt I'm solving this problem in a slightly bad way because I want to
  // mimic (and so test) the code in BrickMeshBemBuilder::CalcFbpsWts. If I ever
  // start to use this routine seriously, I'll optimize the math.
  Matrix<double> A(n, 3), R(3, 3);
  { const double* bx = &_x(1, id+1);
    Matrix<double> Q(n, 3);
    for (size_t j = 0; j < n; j++) {
      RectId src = sids[j];
      Rect sr = _rmu.GetRects()[src];
      if (_any_periodic) {
        Dir::Enum dir;
        double dx, dy;
        Distance(id, src, &dx, &dy, &dir);
        int Dx, Dy;
        Dir::ToVector(dir, &Dx, &Dy);
        sr.x += Dx * _r.dx;
        sr.y += Dy * _r.dy;
      }
      A(j+1, 1) = 1;
      A(j+1, 2) = sr.x - bx[0] + 0.5*sr.dx;
      A(j+1, 3) = sr.y - bx[1] + 0.5*sr.dy;
    }
    SkinnyQr(A, Q, R); }

  double alpha[3];
  { alpha[0] = 1; alpha[1] = 0; alpha[2] = 0;
    blas_int info;
    trtrs('u', 't', 'n', 3, 1, R.GetPtr(), 3, alpha, 3, info);
    assert (info == 0);
    trtrs('u', 'n', 'n', 3, 1, R.GetPtr(), 3, alpha, 3, info);
    assert (info == 0); }
  
  const double* pA = A.GetPtr();
  double z = 0;
  for (size_t j = 0; j < n; j++)
    z += ((alpha[0]*pA[j] + alpha[1]*pA[n + j] + alpha[2]*pA[2*n + j]) *
          z_int[sids[j]]);

#if 0
  // Compare with BrickMeshBemBuilder::CalcFbpsWts debug output.
  Matrix<double> w(n);
  for (size_t j = 0; j < n; j++)
    w(j+1) = alpha[0]*pA[j] + alpha[1]*pA[n + j] + alpha[2]*pA[2*n + j];
  cout << "rws{" << id+1 << "} = " << w << ";" << endl;
#endif

  return z;
}

// -----------------------------------------------------------------------------
// Refine a mesh so that adjacent rectangles differ in size by no more than a
// factor of 2.

class RmuSmoother {
public:
  RmuSmoother (const RectMeshUnstruct* rmu, const RmuAnalyzer* ra)
    : _ormu(rmu), _ora(ra), _srmu(NULL), _rspbl_srmu(true)
  { MakeSmoothRectMeshUnstruct(); }
  ~RmuSmoother () { if (_rspbl_srmu && _srmu) delete _srmu; }

  // User takes responsibility of the pointer after calling this method.
  RectMeshUnstruct* GetSmoothRectMeshUnstruct () {
    _rspbl_srmu = false;
    return _srmu;
  }

private:
  const RectMeshUnstruct* _ormu;
  const RmuAnalyzer* _ora;
  RectMeshUnstruct* _srmu;
  bool _rspbl_srmu;

  struct Node {
    // Each node points to its QuadTree node.
    RectMeshUnstruct::QuadTree* qt;
    // Adjacent rectangles.
    list<Node*> nbrs;
    bool valid;
    Node () : qt(NULL), valid(true) {}
  };
  //todo We might turn this into a vector and remove the lit field from Node.
  list<Node> _nodes; // push_front(), erase()

private:
  void MakeSmoothRectMeshUnstruct();
  void MakeGraph();
  void Smooth();
  void ResetIds();
  bool IsCorrect() const;
  void MakeNodes(RectMeshUnstruct::QuadTree* qt, vector<Node*>& pnodes);
  void Split(Node* node, Node* split_nbrs[nkids]);
  void RemoveNodeFromNbrs(Node* node);
  void SetNbrs(Node* node, list<Node*>& parent_nbrs);
  bool AreAdjacent(const PointRectT<int>& pr, const Rect& nr) const;
  void ResetIds(RectMeshUnstruct::QuadTree* qt, RectId* id);
};

#define foreach(iterator, it, container)                                \
  for (iterator (it) = (container).begin(); (it) != (container).end(); ++(it))

void RmuSmoother::MakeSmoothRectMeshUnstruct () {
  _srmu = new RectMeshUnstruct(*_ormu);
  MakeGraph();
  Smooth();
  ResetIds();
#ifdef AMB_SPECIAL
  // This is an expensive assertion: about 3x as expensive as running the
  // smoother.
  assert(IsCorrect());
#endif
}

void RmuSmoother::MakeGraph () {
  size_t nr = _ormu->GetRects().size();
  vector<Node*> pnodes(nr, NULL);
  for (size_t i = 0; i < _srmu->_qts.size(); i++)
    MakeNodes(_srmu->_qts[i], pnodes);
  for (size_t i = 0; i < nr; i++) {
    Node* node = pnodes[i];
    assert(node);
    assert(node->qt->_is_leaf);
    const vector<RectId>& nbrs = _ora->GetNbrs(node->qt->_l.id);
    for (size_t j = 0; j < nbrs.size(); j++)
      if (nbrs[j] < nr)
        node->nbrs.push_back(pnodes[nbrs[j]]);
  }
#ifndef NDEBUG
  for (size_t i = 0; i < nr; i++) {
    Node* node = pnodes[i];
    foreach (list<Node*>::iterator, nit, node->nbrs) {
      bool fnd = false;
      foreach (list<Node*>::iterator, nnit, (*nit)->nbrs)
        if (*nnit == node) {
          fnd = true;
          break;
        }
      assert(fnd);
    }
  }
#endif
}

void RmuSmoother::
MakeNodes (RectMeshUnstruct::QuadTree* qt, vector<Node*>& pnodes) {
  if (qt->_is_leaf) {
    _nodes.push_front(Node());
    list<Node>::iterator lit = _nodes.begin();
    lit->qt = qt;
    pnodes[qt->_l.id] = &*lit;
  } else {
    for (size_t i = 0; i < nkids; i++)
      MakeNodes(qt->_n.kids[i], pnodes);
  }
}

void RmuSmoother::Smooth () {
  queue<Node*> q;
  // Init q with all rectangles.
  foreach (list<Node>::iterator, it, _nodes) q.push(&*it);
  // Keep checking rectangles until no further splitting is necessary.
  while (!q.empty()) {
    Node* node = q.front();
    q.pop();
    if (!node->valid) continue;
    // Is this rectangle too big?
    bool too_big = false;
    int unit = _ora->UnitDx(node->qt->_r.dx);
    foreach (list<Node*>::iterator, it, node->nbrs) {
      Node* nbr = *it;
      int nunit = _ora->UnitDx(nbr->qt->_r.dx);
      if (unit > 2*nunit) { too_big = true; break; }
    }
    if (too_big) { // Yes.
      // Check the nbr nodes. A node may occupy multiple slots of q; that's ok.
      foreach (list<Node*>::iterator, it, node->nbrs) q.push(*it);
      Node* split_nbrs[nkids];
      Split(node, split_nbrs);
      // Check the new nodes.
      for (size_t i = 0; i < nkids; i++) q.push(split_nbrs[i]);
      // We can't erase the list entry because q might still have it, so just
      // invalidate it.
      node->valid = false;
    }
  }
}

void RmuSmoother::Split (Node* node, Node* split_nbrs[nkids]) {
  RemoveNodeFromNbrs(node);
  node->qt->Split();
  for (size_t i = 0; i < nkids; i++) {
    _nodes.push_front(Node());
    list<Node>::iterator lit = _nodes.begin();
    lit->qt = node->qt->_n.kids[i];
    SetNbrs(&*lit, node->nbrs);
    split_nbrs[i] = &*lit;
  }
  for (size_t i = 0; i < nkids; i++)
    for (size_t j = 0; j < nkids; j++) {
      if (i == j) continue;
      split_nbrs[i]->nbrs.push_back(split_nbrs[j]);
    }
}

void RmuSmoother::RemoveNodeFromNbrs (Node* node) {
  foreach (list<Node*>::iterator, nit, node->nbrs) {
    Node* nbr = *nit;
    bool fnd = false;
    foreach (list<Node*>::iterator, nnit, nbr->nbrs)
      if (*nnit == node) {
        nbr->nbrs.erase(nnit);
        fnd = true;
        break;
      }
    assert(fnd);
  }
}

void RmuSmoother::SetNbrs (Node* node, list<Node*>& pnbrs) {
  const Rect& r = node->qt->_r;
  PointRectT<int> pr(RectT<int>(_ora->UnitX(r.x), _ora->UnitY(r.y),
                                _ora->UnitDx(r.dx), _ora->UnitDy(r.dy)));
  foreach (list<Node*>::iterator, pnit, pnbrs) {
    Node* pnbr = *pnit;
    if (AreAdjacent(pr, pnbr->qt->_r)) {
      node->nbrs.push_back(pnbr);
      pnbr->nbrs.push_back(node);
    }
  }
}

bool RmuSmoother::
AreAdjacent (const PointRectT<int>& pr, const Rect& nr) const {
  PointRectT<int> npr(RectT<int>(_ora->UnitX(nr.x), _ora->UnitY(nr.y),
                                 _ora->UnitDx(nr.dx), _ora->UnitDy(nr.dy)));
  // Easy exit: adjacent as is.
  if (PointRectT<int>::IntersectClosed(pr, npr)) return true;
  bool per_N = _ora->GetBoundaries().GetBC(Dir::N) == Boundaries::bc_periodic;
  bool per_E = _ora->GetBoundaries().GetBC(Dir::E) == Boundaries::bc_periodic;

  // Another easy exit: no periodic BCs to worry about.
  if (!(per_N || per_E)) return false;

  // Now we need to move nr. Either there are 2 places to check or 8.
  RectT<int> nr_unit(_ora->UnitX(nr.x), _ora->UnitY(nr.y),
                     _ora->UnitDx(nr.dx), _ora->UnitDy(nr.dy));
#define check(delta_E, delta_N) do {                                    \
    { PointRectT<int> npr(nr_unit, delta_E, delta_N);                   \
      if (PointRectT<int>::IntersectClosed(pr, npr)) return true; }     \
  } while (0)
  int delta_E = _ora->UnitDx(_ormu->GetDomain().dx);
  if (per_E) {
    check( delta_E, 0);
    check(-delta_E, 0);
  }
  int delta_N = _ora->UnitDy(_ormu->GetDomain().dy);
  if (per_N) {
    check(0,  delta_N);
    check(0, -delta_N);
    if (per_E) {
      check( delta_E,  delta_N);
      check( delta_E, -delta_N);
      check(-delta_E,  delta_N);
      check(-delta_E, -delta_N);
    }
  }
  return false;
#undef check
}

void RmuSmoother::ResetIds () {
  RectId id = 0;
  for (size_t i = 0; i < _srmu->_qts.size(); i++)
    ResetIds(_srmu->_qts[i], &id);
}

void RmuSmoother::ResetIds (RectMeshUnstruct::QuadTree* qt, RectId* id) {
  if (qt->_is_leaf) qt->_l.id = (*id)++;
  else for (size_t i = 0; i < nkids; i++) ResetIds(qt->_n.kids[i], id);
}

template<typename T>
ostream& operator<< (ostream& os, const vector<T>& v) {
  os << "[ ";
  for (size_t i = 0; i < v.size(); i++) os << v[i] << " ";
  return os << "]";
}

bool RmuSmoother::IsCorrect () const {
  const vector<Rect>& rs = _srmu->GetRects();
  size_t nr = rs.size();
  RmuAnalyzer ra(_srmu, _ora->GetBoundaries());
  foreach (list<Node>::const_iterator, it, _nodes) {
    if (!it->valid) continue;
    if (!it->qt->_is_leaf) {
      cerr << "Not a leaf." << endl;
      return false;
    }
    int unit = _ora->UnitDx(it->qt->_r.dx);
    foreach (list<Node*>::const_iterator, nit, it->nbrs) {
      int nunit = _ora->UnitDx((*nit)->qt->_r.dx);
      if (unit > 2*nunit || nunit > 2*unit) {
        cerr << "Not smooth." << endl;
        return false;
      }
    }
    const vector<RectId>& nbrs = ra.GetNbrs(it->qt->_l.id);
    vector<RectId> nids1, nids2;
    for (size_t i = 0; i < nbrs.size(); i++)
      if (nbrs[i] < nr) nids1.push_back(nbrs[i]);
    foreach (list<Node*>::const_iterator, nit, it->nbrs)
      nids2.push_back((*nit)->qt->_l.id);
    sort(nids1.begin(), nids1.end());
    sort(nids2.begin(), nids2.end());
    if (nids1.size() != nids2.size()) {
      fprintf(stderr, "|nids1| %ld |nids2| %ld\n", nids1.size(), nids2.size());
      cerr << "nids1: " << nids1 << endl;
      cerr << "nids2: " << nids2 << endl;
      cout << "master: " << rs[it->qt->_l.id] << endl;
      for (size_t i = 0; i < nids2.size(); i++)
        cout << rs[nids2[i]] << endl;
      return false;
    }
    for (size_t i = 0; i < nids1.size(); i++)
      if (nids1[i] != nids2[i]) {
        cerr << "Different nbr ids." << endl;
        return false;
      }
  }
  return true;
}

#undef foreach

RectMeshUnstruct*
SmoothRectMeshUnstruct (const RectMeshUnstruct* rmu, const Boundaries& b) {
  RmuAnalyzer* ra = new RmuAnalyzer(rmu, b, true); // Adjacency graph only.
  RmuSmoother* rs = new RmuSmoother(rmu, ra);
  RectMeshUnstruct* srmu = rs->GetSmoothRectMeshUnstruct();
  delete rs;
  DeleteRmuAnalyzer(ra);
  return srmu;
}
}}
