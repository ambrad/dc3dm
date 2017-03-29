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
namespace es {

inline double GetTractionComp (
  const LameParms& lp, const dc3::Elem& es, const double* disl,
  const dc3::Elem& eo, size_t component, bool want_fullspace)
{
  double u[3], du[9], s[6];
  Dc3d(lp, es, disl, eo.Center(), u, du, want_fullspace);
  DuToS(lp, du, s);
  return ProjectStress(s, eo, component);
}

inline double
ProjectStress (const double* s, const dc3::Elem& e, size_t component) {
  double B;
  switch (component) {
  case 0:
    ProjectStress(s, e.Normal(), e.AlongStrike(), NULL, &B,   NULL, NULL);
    break;
  case 1:
    ProjectStress(s, e.Normal(), e.AlongDip(),    NULL, &B,   NULL, NULL);
    break;
  case 2:
    ProjectStress(s, e.Normal(), NULL,            NULL, NULL, NULL, &B  );
    break;
  }
  return B;      
}

inline void dc3::Elem::ToGlobal (const double vecl[3], double vecg[3]) const {
  vecg[0] = _astrike[0]*vecl[0] + _adip[0]*vecl[1] + _normal[0]*vecl[2];
  vecg[1] = _astrike[1]*vecl[0] + _adip[1]*vecl[1] + _normal[1]*vecl[2];
  vecg[2] = _astrike[2]*vecl[0] + _adip[2]*vecl[1] + _normal[2]*vecl[2];
}
    
}}
