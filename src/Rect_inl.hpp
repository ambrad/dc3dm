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
 
#ifndef INCLUDE_UTIL_RECT_INL
#define INCLUDE_UTIL_RECT_INL

#include "util/include/CodeAnalysis.hpp"
#include "util/include/Util.hpp"

namespace util {
namespace rmesh {
inline Dir::Enum Dir::ToCardinal (Dir::Enum dir) {
  switch (dir) {
  case Dir::SE: case Dir::NE: return Dir::E;
  case Dir::SW: case Dir::NW: return Dir::W;
  default: return dir;
  }
}

inline bool Dir::IsCardinal (Dir::Enum dir) { return (int) dir < 4; }

inline bool Dir::IsCorner(Dir::Enum dir) { return (int) dir >= 5; }

inline void Dir::BreakUp (Dir::Enum dir, Dir::Enum& dir1, Dir::Enum& dir2) {
  dir1 = dir2 = inside;
  if (IsCardinal(dir)) {
    switch (dir) {
    case N: case S: dir1 = dir; break;
    case E: case W: dir2 = dir; break;
    default: break;
    }
    return;
  }
  switch (dir) {
  case NE: dir1 = N; dir2 = E; break;
  case NW: dir1 = N; dir2 = W; break;
  case SW: dir1 = S; dir2 = W; break;
  case SE: dir1 = S; dir2 = E; break;
  default: assert(false);
  }
}

inline Dir::Enum Dir::Combine (Dir::Enum dir1, Dir::Enum dir2) {
  assert((int) dir1 < 5 && (int) dir2 < 5); // cardinal or inside
  assertpr(dir1 == inside || dir1 != Opposite(dir2),
           "dir1 = %d  dir2 = %d", (int) dir1, (int) dir2);
  if ((int) dir2 < (int) dir1) std::swap(dir1, dir2);
  if (dir1 == E) {
    if (dir2 == N) return NE;
    if (dir2 == S) return SE;
  }
  if (dir1 == N && dir2 == W) return NW;
  if (dir1 == W && dir2 == S) return SW;
  assert(dir2 == Dir::inside || dir2 == dir1);
  assert(dir1 == Dir::inside || IsCardinal(dir1));
  return dir1;
}

inline Dir::Enum Dir::Opposite (Dir::Enum side) {
  switch (side) {
  case E:  return W;
  case NE: return SW;
  case N:  return S;
  case NW: return SE;
  case W:  return E;
  case SW: return NE;
  case S:  return N;
  case SE: return NW;
  default: return inside;
  }
}

inline Dir::Enum Dir::OneCw (Dir::Enum dir) {
  switch (dir) {
  case E:  return SE;
  case NE: return E;
  case N:  return NE;
  case NW: return N;
  case W:  return NW;
  case SW: return W;
  case S:  return SW;
  case SE: return S;
  default: return inside;
  }
}

inline Dir::Enum Dir::OneCcw (Dir::Enum dir) {
  switch (dir) {
  case E:  return NE;
  case NE: return N;
  case N:  return NW;
  case NW: return W;
  case W:  return SW;
  case SW: return S;
  case S:  return SE;
  case SE: return E;
  default: return inside;
  }
}

inline Dir::Enum Dir::ToEnum (int Dx, int Dy) {
  Dir::Enum dir1, dir2;
  switch (Dy) {
  case -1: dir1 = S; break;
  case  0: dir1 = inside; break;
  case  1: dir1 = N; break;
  default: assertpr(false, "Dy = %d", Dy);
  }
  switch (Dx) {
  case -1: dir2 = W; break;
  case  0: dir2 = inside; break;
  case  1: dir2 = E; break;
  default: assertpr(false, "Dx = %d", Dx);
  }
  return Combine(dir1, dir2);
}

inline void Dir::ToVector (Dir::Enum dir, int* Dx, int* Dy) {
  switch (dir) {
  case E:  *Dx =  1; *Dy =  0; break;
  case NE: *Dx =  1; *Dy =  1; break;
  case N:  *Dx =  0; *Dy =  1; break;
  case NW: *Dx = -1; *Dy =  1; break;
  case W:  *Dx = -1; *Dy =  0; break;
  case SW: *Dx = -1; *Dy = -1; break;
  case S:  *Dx =  0; *Dy = -1; break;
  case SE: *Dx =  1; *Dy = -1; break;
  default: *Dx = 0; *Dy = 0; break;
  }
}

inline Boundaries::Boundaries () {
  for (size_t i = 0; i < 4; i++) _bcs[i] = bc_0velocity;
}

inline void Boundaries::Serialize (FILE* fid) const throw (FileException) {
  write(_bcs, 4, fid);
  write(_rs, 4, fid);
}

inline void Boundaries::Deserialize (FILE* fid) throw (FileException) {
  read(_bcs, 4, fid);
  read(_rs, 4, fid);
}

inline void Boundaries::SetVelocityBC (Dir::Enum side, const Rect& r) {
  side = Dir::ToCardinal(side);
  _bcs[side] = bc_velocity;
  _rs[side] = r;
}

inline void Boundaries::SetZeroVelocityBC (Dir::Enum side) {
  _bcs[Dir::ToCardinal(side)] = bc_0velocity;
}

inline void Boundaries::SetPeriodicBC (Dir::Enum side) {
  side = Dir::ToCardinal(side);
  _bcs[side] = bc_periodic;
  _bcs[Dir::Opposite(side)] = bc_periodic;
}

inline void Boundaries::SetFreeBC (Dir::Enum side) throw (Exception) {
  size_t cnt = 0;
  for (size_t i = 0; i < 4; i++) if (_bcs[i] == bc_free) cnt++;
  if (cnt == 1)
    throw Exception("Boundaries: One BC is already set to bc_free.");
  assert(cnt == 0);
  _bcs[side] = bc_free;
}

inline Boundaries::BC Boundaries::GetBC (Dir::Enum side) const {
  int idx = (int) Dir::ToCardinal(side);
  assert(idx >= 0 && idx < 4);
  return _bcs[idx];
}

inline const Rect& Boundaries::GetRect (Dir::Enum side) const {
  int idx = (int) Dir::ToCardinal(side);
  assert(idx >= 0 && idx < 4);
  return _rs[idx];
}

inline bool Boundaries::IsVBC (Boundaries::BC bc) {
  return bc == bc_velocity || bc == bc_0velocity || bc == bc_free;
}

inline bool Boundaries::HasVBC (Dir::Enum dir) const {
  if (dir == Dir::inside) return false;
  if (Dir::IsCardinal(dir)) return IsVBC(GetBC(dir));
  else {
    Dir::Enum dir1, dir2;
    Dir::BreakUp(dir, dir1, dir2);
    return IsVBC(GetBC(dir1)) || IsVBC(GetBC(dir2));
  }
}

inline bool Boundaries:: HasZeroOrOneFree() const {
  size_t cnt = 0;
  for (size_t i = 0; i < 4; i++) if (_bcs[i] == bc_free) cnt++;
  return cnt <= 1;
}
}}

#endif
