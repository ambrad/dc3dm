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
 
namespace util {
namespace rmesh {
using namespace std;

// NB: Unit coords here are half those in RmuAnalyzer.
inline int MeshAnalyzer::UnitDx (double dx) const {
  int idx = (int) round(dx / _ux);
  assertpr(fabs(idx - dx / _ux) / (1 + fabs(idx)) <= 1e-10,
           " _ux = %1.2f dx = %1.2f idx = %d", _ux, dx, idx);
  return idx;
}

inline int MeshAnalyzer::UnitDy (double dy) const {
  int idy = (int) round(dy / _uy);
  assertpr(fabs(idy - dy / _uy) / (1 + fabs(idy)) <= 1e-10,
           " _uy = %1.2f dy = %1.2f idy = %d", _uy, dy, idy);
  return idy;
}

template<typename T> inline T Square (T x) { return x*x; }

inline double MeshAnalyzer::
Distance (RectId id1, RectId id2, Dir::Enum* dir) const {
  double dx, dy;
  Distance(id1, id2, &dx, &dy, dir);
  assert (GetBoundaries().GetBC(Dir::E) != Boundaries::bc_periodic
          || UnitDx(dx) < 0.5 * UnitDx(GetDomain().dx));
  assert (GetBoundaries().GetBC(Dir::N) != Boundaries::bc_periodic
          || UnitDy(dy) < 0.5 * UnitDy(GetDomain().dy));
#ifdef EUCLIDEAN_DISTANCE
  return sqrt((double) Square(UnitDx(dx)) + (double) Square(UnitDy(dy)));
#else
  return (size_t) std::max(UnitDx(dx), UnitDy(dy));
#endif
}

inline double MeshAnalyzer::
Distance (const Rect& r1, const Rect& r2) const {
  double dx, dy;
  Rect::Distance(r1, r2, &dx, &dy);
#ifdef EUCLIDEAN_DISTANCE
  return sqrt((double) Square(UnitDx(dx)) + (double) Square(UnitDy(dy)));
#else
  return (size_t) std::max(UnitDx(dx), UnitDy(dy));
#endif
}

inline Rect& MeshAnalyzer::TranslateRect (Rect& r, Dir::Enum dir) const {
  if (dir == Dir::inside) return r;
  int Dx, Dy;
  Dir::ToVector(dir, &Dx, &Dy);
  return TranslateRect(r, Dx, Dy);
}

inline Rect& MeshAnalyzer::TranslateRect (Rect& r, int Dx, int Dy) const {
  r.x += Dx * GetDomain().dx;
  r.y += Dy * GetDomain().dy;
  return r;
}

}}
