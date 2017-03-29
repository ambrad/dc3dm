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
 
#ifndef INCLUDE_UTIL_RECT
#define INCLUDE_UTIL_RECT

#include "util/include/Exception.hpp"

namespace util {
namespace rmesh {

// ---------------------------------------------------------------------------
// Describe a rectangle domain and its boundary conditions.

typedef size_t RectId;

template<typename T>
struct RectT {
  T x, y;   // lower-left corner
  T dx, dy; // x,y-direction lengths

  RectT() : x(0.0), y(0.0), dx(0.0), dy(0.0) {}
  RectT(T ix, T iy, T idx, T idy) : x(ix), y(iy), dx(idx), dy(idy) {}

  static void Distance(const RectT<T>& r1, const RectT<T>& r2, T* dx, T* dy);
};

typedef RectT<double> Rect;

struct Dir {
  enum Enum { E = 0, N, W, S, inside, NE, NW, SW, SE };

  // Implements conventions for convering 8-compass-direction dir to cardinal
  // directions plus inside.
  static Enum ToCardinal(Enum dir);
  static bool IsCardinal(Enum dir);
  static bool IsCorner(Enum dir);
  // Convert between non-cardinal and 2 cardinal directions, if necessary. Pad
  // with Dir::inside.
  //   On exit, dir1 is N, S, or inside; dir2 is E, W, or inside.
  static void BreakUp(Enum dir, Enum& dir1, Enum& dir2);
  //   On entry, dir1, dir2 can each be cardinal or inside.
  static Enum Combine(Enum dir1, Enum dir2);
  // Encode a direction as the integer-valued vector (\pm Dx, \pm Dy).
  static Enum ToEnum(int Dx, int Dy);
  static void ToVector(Enum dir, int* Dx, int* Dy);
  // Opposite direction.
  static Enum Opposite(Enum dir);
  // Go one clockwise and...
  static Enum OneCw(Enum dir);
  // ...counter-clockwise. For example: E -> NE.
  static Enum OneCcw(Enum dir);
};

// Set the east-west and north-south boundary conditions. All directions are
// converted to cardinal.
class Boundaries {
public:
  enum BC { bc_velocity, bc_0velocity, bc_periodic, bc_free };

  // By default, every boundary has boundary condition bc_0velocity.
  Boundaries();

  void Serialize(FILE* fid) const throw (FileException);
  void Deserialize(FILE* fid) throw (FileException);

  void SetVelocityBC
  (Dir::Enum side,
   // A velocity BC is implemented by a large rectangle. To be sensible, the
   // rectangle should be adjacent to the side.
   const Rect& r);
  // A 0-velocity BC can be treated specially for efficiency.
  void SetZeroVelocityBC(Dir::Enum side);
  // The opposite side is also set since that is necessary for periodic BC.
  void SetPeriodicBC(Dir::Enum side);
  // Only one side can be bc_free. If another side is already set to bc_free,
  // an exception is thrown.
  void SetFreeBC(Dir::Enum side) throw (Exception);

  BC GetBC(Dir::Enum side) const;
  const Rect& GetRect(Dir::Enum side) const;

  // From the perspective of mesh operations, all BCs except bc_periodic act
  // the same.
  static bool IsVBC(BC bc);
  // If dir is a corner, return true if either side of the corner is a v-BC.
  bool HasVBC(Dir::Enum dir) const;

  bool HasZeroOrOneFree() const;

private:
  BC _bcs[4]; // ordered _bcs[Dir::Enum]
  Rect _rs[4];
};

}}

#endif
#include "Rect_inl.hpp"
