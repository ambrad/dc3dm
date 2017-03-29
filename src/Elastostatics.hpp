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
 
#ifndef INCLUDE_DC3DM_ELASTOSTATICS
#define INCLUDE_DC3DM_ELASTOSTATICS

#include <math.h>
#include <vector>
#include "util/include/Matrix.hpp"
#include "util/include/Exception.hpp"
#include "RectMeshUnstruct_pri.hpp"

namespace util {
using namespace std;
using namespace util;

namespace es {

inline double cosd (double a) { return cos(a*M_PI/180); }
inline double sind (double a) { return sin(a*M_PI/180); }

class LameParms {
public:
  LameParms () {}
  LameParms (double mu, double nu) { Set(mu, nu); }

  void Set (double mu,   // Shear modulus
            double nu) { // Poisson's ratio
    _mu = mu;
    _nu = nu;
    _lambda = 2*mu*nu/(1 - 2*nu);
    _alpha  = (_lambda + mu)/(_lambda + 2.0*mu);
  }

  double mu ()     const { return _mu;     }
  double nu ()     const { return _nu;     }
  double lambda () const { return _lambda; }
  double alpha ()  const { return _alpha;  }
      
private:
  double _mu, _nu, _lambda, _alpha;
};

namespace dc3 {
// Various interfaces to Y. Okada's dc3d routine.

// Declaration for the Fortran routine. The rectangle slants upward for
// positive dip in the y direction.
extern "C" void dc3d_(
  const char* space, const double* ALPHA, const double* X, const double* Y,
  const double* Z, const double* DEPTH, const double* DIP, const double* AL1,
  const double* AL2, const double* AW1, const double* AW2,
  const double* DISL1, const double* DISL2, const double* DISL3,
  double* UX, double* UY, double* UZ, double* UXX, double* UYX,
  double* UZX, double* UXY, double* UYY, double* UZY, double* UXZ,
  double* UYZ, double* UZZ);

inline void Dc3d (
  double alpha,
  double x, double y, double z, // NB: x, y are relative to the rectangle
  double depth, double dipdeg,
  double al1, double al2, double aw1, double aw2,
  double disl_strike, double disl_dip, double disl_tensile,
  double* u,  // [ux uy uz]
  double* du, // [uxx uyx uzx uxy uyy uzy uxz uyz uzz]
  bool want_fullspace = false)
{
  char wsh = 'H';
  if (want_fullspace) wsh = 'f';
  dc3d_(&wsh, &alpha,
        &x, &y, &z,
        &depth, &dipdeg, &al1, &al2, &aw1, &aw2,
        &disl_strike, &disl_dip, &disl_tensile,
        u, u+1, u+2,
        du, du+1, du+2, du+3, du+4, du+5, du+6, du+7, du+8);
}

class Elem {
public:
  Elem(double depth, double dipdeg,
       double al1, double al2,
       double aw1, double aw2,
       double gx, double gy);
  Elem(double depth,
       double dipdeg, double cos_dip, double sin_dip, // avoid trig
       double al1, double al2,
       double aw1, double aw2,
       double gx, double gy);
  void Set(double depth, double dipdeg, double cos_dip, double sin_dip,
           double al1, double al2, double aw1, double aw2,
           double gx, double gy);

  // dc3d element specification:
  double depth ()  const { return _depth;  }
  double dipdeg () const { return _dipdeg; }
  // x, along-strike
  double al1 () const { return _al1; }
  double al2 () const { return _al2; }
  // y, along-dip
  double aw1 () const { return _aw1; }
  double aw2 () const { return _aw2; }

  // Global location of element's local origin. Use these to transform a
  // global obs vector to this element's coordinate system before calling
  // dc3d.
  double gx () const { return _gx; }
  double gy () const { return _gy; }

  // Global center of this element [x y z] (NB: z, not depth):
  const double* Center ()      const { return _ctr;     }
  // Unit in-plane vector pointing in the y direction
  const double* AlongDip ()    const { return _adip;    }
  // Unit in-plane vector pointing in the x direction
  const double* AlongStrike () const { return _astrike; }
  // Unit element normal vector
  const double* Normal ()      const { return _normal;  }

  // Transform from a vector specifed in the [strike, dip, normal] coordinate
  // system to global [x, y, z] coordinates.
  void ToGlobal(const double vec_local[3], double vec_global[3]) const;

private:
  double _depth, _dipdeg, _al1, _al2, _aw1, _aw2;
  double _gx, _gy;
  double _ctr[3], _adip[3], _astrike[3], _normal[3];
};

inline void Dc3d (
  const LameParms& lp,
  const Elem& e,
  const double* disl, // [strike dip tensile]
  const double* obs,  // [x y z] NB: global coordinates
  double* u,          // [ux uy uz]
  double* du,         // [uxx uyx uzx uxy uyy uzy uxz uyz uzz]
  bool want_fullspace = false) {
  Dc3d(lp.alpha(),
       obs[0] - e.gx(), obs[1] - e.gy(), obs[2],
       e.depth(), e.dipdeg(), e.al1(), e.al2(), e.aw1(), e.aw2(),
       disl[0], disl[1], disl[2],
       u, du, want_fullspace);
}

// Create a vector of elements, ordered with fastest index -x to +x, from
// a tensor-mesh specification of a planar, rectangular fault.
void PlanarTensorMeshToElems(
  const Matrix<double>& x,   // Along-strike element edges in increasing x
  const Matrix<double>& eta, // Along-dip, in-plane edges in increasing eta
  double y_min,              // min y value of mesh
  double depth_min,          // min depth of mesh
  double dipdeg,             // Dip angle
  vector<Elem>& es)          // NB: es is *not* cleared prior to filling it
  throw (Exception); // If the mesh sticks out of the ground

void RectMeshUnstructToElems(
  const rmesh::RectMeshUnstruct& rm,
  double y_min, double depth_min, double dipdeg, vector<Elem>& es)
  throw (Exception);
} // namespace dc3
        
void DuToS(const LameParms& lp,
           const double* du, // du as in dc3
           double* s);       // [sxx sxy sxz syy syz szz]

// From s, compute the shear and normal stresses. This routine does not
// check for math errors.
void ProjectStress(
  const double* s,      // 6-element output from DuToS that gives Sigma.
  const double* normal, // 3-element unit vector normal to the element.
  const double* along1, // Optional 3-element unit vector along the element
  const double* along2, // and optional second one; set to NULL otherwise.
  double* ss1,          // along1'*Sigma*normal if along1 is not NULL
  double* ss2,          // and along2'*Sigma*normal if along2 is not NULL.
  double* ns);          // Optionally request normal(:,i)'*Sigma*normal.

// Convenience wrapper to Dc3d, DuToS, ProjectStress. Return one component
// of traction due to the influence of the source element es on the center
// of the observation element eo.
double GetTractionComp(
  const LameParms& lp,
  const dc3::Elem& es, // source
  const double* disl,  // [strike dip tensile]
  const dc3::Elem& eo, // eo.Center() is obs
  size_t component,    // 0 - strike, 1 - dip, 2 - normal
  bool want_fullspace = false);

double ProjectStress(const double* s, const dc3::Elem& e, size_t component);
} // namespace es
} // namespace util

#include "Elastostatics_inl.hpp"

#endif
