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
 
#ifndef INCLUDE_UTIL_TRI2
#define INCLUDE_UTIL_TRI2

#include <math.h>
#include <assert.h>

namespace util {
namespace rmesh {

// Low-level routiens for dealing with triangles.

class Tri2 {
  double const* _tri[3];

public:
  Tri2 (const double a[2], const double b[2], const double c[2])
  { Set(a, b, c); }
  // Use at your own risk.
  Tri2 () { Set(NULL, NULL, NULL); }

  void Set (const double a[2], const double b[2], const double c[2])
  { _tri[0] = a; _tri[1] = b; _tri[2] = c; }

  static size_t Nvi (size_t i) { return (i + 1) % 3; }
  static size_t Pvi (size_t i) { return (i + 2) % 3; }

  static bool IsCcw (const double* a, const double* b, const double* c) {
    // Is cross([b - a; 0], [c - a; 0]) >= 0?
    return (b[0] - a[0])*(c[1] - a[1]) >= (b[1] - a[1])*(c[0] - a[0]);
  }
  bool IsCcw () const { return IsCcw(_tri[0], _tri[1], _tri[2]); }

  double EdgeLen2 (size_t i, size_t j) const {
    return Square(_tri[j][0] - _tri[i][0]) +
      Square(_tri[j][1] - _tri[i][1]);
  }
  double MaxEdgeLen2(size_t* max_edge_len_idx = NULL) const;

  // Triangle's aspect ratio.
  double AspectRatio(size_t* max_edge_len_idx = NULL) const;

  void BarycentricMatrix(double Ti[4]) const;
  void ToBarycentric(const double Ti[4], const double x[2], double lambda[3])
    const;
  double ToBarycentric(const double Ti[4], const double x[2],
                       // Get just this componet of lambda.
                       size_t i_lambda) const;

  template<typename T> static inline T Square (T v) { return v*v; }
  // c = a - b
  static double* Diff (const double a[2], const double b[2], double c[2])
  { c[0] = a[0] - b[0]; c[1] = a[1] - b[1]; return c; }
  void Diff (size_t i, size_t j, double d[2]) const
  { Diff(_tri[i], _tri[j], d); }
  static double* Scale (double v[2], double a)
  { v[0] *= a; v[1] *= a; return v; }
  static double Dot (const double a[2], const double b[2])
  { return a[0]*b[0] + a[1]*b[1]; }
  static double Norm2 (const double v[2]) { return Dot(v, v); }
};

inline double Tri2::MaxEdgeLen2 (size_t* max_edge_len_idx) const {
  size_t meli;
  double mel = 0.0;
  for (int i = 0; i < 3; i++) {
    double len2 = EdgeLen2(i, Nvi(i));
    if (len2 > mel) { mel = len2; meli = i; }
  }
  if (max_edge_len_idx) *max_edge_len_idx = meli;
  return mel;
}

inline double Tri2::AspectRatio (size_t* max_edge_len_idx) const {
  size_t meli;
  double mel = MaxEdgeLen2(&meli);

  double perp;
  { double b[2], v[2];
    Diff(Nvi(meli), meli, b);
    Diff(Pvi(meli), meli, v);
    Diff(v, Scale(b, Dot(b, v)/mel), v);
    perp = Norm2(v); }

  if (max_edge_len_idx) *max_edge_len_idx = meli;
  return sqrt(mel / perp);
}

inline void Tri2::BarycentricMatrix (double Ti[4]) const {
  Diff(_tri[0], _tri[2], Ti);
  Diff(_tri[1], _tri[2], Ti + 2);
  const double det = Ti[0]*Ti[3] - Ti[1]*Ti[2];
  const double tmp = Ti[0];
  Ti[0] = Ti[3] / det;
  Ti[1] = -Ti[1] / det;
  Ti[2] = -Ti[2] / det;
  Ti[3] = tmp / det;
}

inline void Tri2::
ToBarycentric (const double Ti[4], const double x[2], double lam[3]) const {
  double dx[2];
  Diff(x, _tri[2], dx);
  lam[0] = Ti[0] * dx[0] + Ti[2] * dx[1];
  lam[1] = Ti[1] * dx[0] + Ti[3] * dx[1];
  lam[2] = 1.0 - lam[0] - lam[1];
}

inline double Tri2::
ToBarycentric (const double Ti[4], const double x[2], size_t i_lambda) const {
  double dx[2];
  Diff(x, _tri[2], dx);
  switch (i_lambda) {
  case 0: return Ti[0] * dx[0] + Ti[2] * dx[1];
  case 1: return Ti[1] * dx[0] + Ti[3] * dx[1];
  case 2:
    return 1.0 - (Ti[0] * dx[0] + Ti[2] * dx[1]) -
      (Ti[1] * dx[0] + Ti[3] * dx[1]);
  default:
    assert(false);
    return 0;
  }
}
  
}}

#endif
