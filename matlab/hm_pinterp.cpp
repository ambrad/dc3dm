/* hmmvp: Software to form and apply Hierarchical Matrices
 *   Version 1.3
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 *   https://pangea.stanford.edu/research/CDFM/software
 * hmmvp is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.0
*/

//$ mex -I../.. hm_pinterp.cpp /home/ambrad/code/util/src/PolyInterp.cpp

#include "MexUtil.hpp"
#include "util/include/PolyInterp.hpp"
using namespace util;
using namespace Interp;

PolyMesh* LToPolyMesh (const mxArray* mL) throw (Exception) {
  Matd xlims;
  MexToMatrix(mxGetField(mL, 0, "xlims"), xlims);
  int nd = xlims.Size(1);
  vector<Matd> x(nd), w(nd);
  const mxArray* Lx = mxGetField(mL, 0, "x");
  const mxArray* Lw = mxGetField(mL, 0, "w");
  for (int i = 0; i < nd; i++) {
    MexToMatrix(mxGetCell(Lx, i), x[i]);
    MexToMatrix(mxGetCell(Lw, i), w[i]);
  }
  return new PolyMesh(xlims, x, w);
}

mxArray* PolyMeshToL (const PolyMesh& pm, const Matd& F) {
  const char* flds[] = {"ns", "xlims", "x", "w", "f"};
  mwSize nL = 1;
  mxArray* mL = mxCreateStructArray(1, &nL, 5, flds);
  const vector<Matd>& x = pm.GetX();
  const vector<Matd>& w = pm.GetW();
  mwSize nd = x.size();
  Vecui ns(nd);
  mxArray* mx = mxCreateCellArray(1, &nd);
  mxArray* mw = mxCreateCellArray(1, &nd);
  for (uint di = 0; di < nd; di++) {
    ns[di] = x[di].Size();
    mxSetCell(mx, di, MatrixToMex(x[di]));
    mxSetCell(mw, di, MatrixToMex(w[di]));
  }
  SetField(mL, 0, flds[0], ns);
  SetField(mL, 0, flds[1], pm.GetXlims());
  mxSetField(mL, 0, flds[2], mx);
  mxSetField(mL, 0, flds[3], mw);
  SetField(mL, 0, flds[4], F);
  return mL;
}

void EvalUserFn (const mxArray* usrfn, const Matd& X, Matd& F) throw (Exception)
{
  int nlhs = 1, nrhs = 2;
  mxArray *prhs[2], *plhs[1], *mX;
  // mexCallMATLAB is supposed to let prhs be const, but it isn't.
  mxArray* usrfn1 = mxDuplicateArray(usrfn);
  prhs[0] = usrfn1;
  mX = MatrixToMex(X);
  prhs[1] = mX;
  mxArray* merr = mexCallMATLABWithTrap(nlhs, plhs, nrhs, prhs, "feval");
  if (merr) throw Exception("Error in user function.");
  mxDestroyArray(mX);
  MexToMatrix(plhs[0], F);
  mxDestroyArray(usrfn1);
}

/* This is a rather messily implemented interface to PolyInterp for my own
   use. The functions are not consistent with each other.

   Implements

       L = hm_pinterp('init',xlims,ns,usrfnstr,[L0])
         xlims is an nd x 2 array of endpoints.
         ns(i) is the number of nodes in dimension i.
         fh is a function handle such that F = fh(X), where X is nxd and each
           row is a d-vector of coordinates at which to evaluate F.
         Note: This function is not really used in practice.

       Fi = hm_pinterp('eval',L,F,Xi)
         Note: This function tends to be used by itself.
         Xi is an n x nd array of points at which to interpolate.
         F replaces L.f. It can have ncomp columns, where ncomp is the number
           of components for F a vector-vauled function.
*/
void mexFunction (int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
  if (nrhs < 1)
    mexErrMsgTxt("First argument must be one of: init, eval.");

  char fn[2];
  mxGetString(prhs[0], fn, 2);
  switch (fn[0]) {
  // init with cheb nodes
  case 'i':
  case 'I': {
    if (nrhs < 4 || nrhs > 5 || nlhs != 1)
      mexErrMsgTxt("L = hm_pinterp('init',xlims,ns,usrfnstr,[L0])");

    Matd xlims;
    MexToMatrix(prhs[1], xlims);
    Vecui ns;
    MexToVector(prhs[2], ns);
    ChebMesh cm(xlims, ns);

    //if (nrhs > 4) mexErrMsgTxt("Not implemented yet.");
    Matd F;
    try {
      Matd X;
      cm.GetNodes(X);
      EvalUserFn(prhs[3], X, F);
    } catch(Exception& e) {
      mexErrMsgTxt(e.GetMsg().c_str());
    }
    plhs[0] = PolyMeshToL(cm, F);

    break;
  }

  // eval
  case 'e':
  case 'E': {
    if (nrhs != 4 || nlhs != 1)
      mexErrMsgTxt("Fi = hm_pinterp('eval',L,F,Xi)");

    PolyMesh* pm = 0;
    try {
      pm = LToPolyMesh(prhs[1]);
    } catch (Exception& e) {
      mexErrMsgTxt(e.GetMsg().c_str());
    }
    Matd F, Xi, Fi;
    MexToMatrix(prhs[2], F);
    MexToMatrix(prhs[3], Xi);

    int nx = Xi.Size(1), nd = Xi.Size(2), ncomp = F.Size(2);
    Fi.Resize(nx,ncomp);
    Matd xi(nd), wrk, Fe(ncomp);
    for (int i = 1; i <= nx; i++) {
      for (int j = 1; j <= nd; j++) xi(j) = Xi(i,j);
#if 1
      // Allow multiple components.
      pm->Eval(F, xi, wrk, Fe);
      for (int j = 1; j <= ncomp; j++) Fi(i,j) = Fe(j);
#else
      Fi(i) = pm->Eval(F, xi, wrk);
#endif
    }
    plhs[0] = MatrixToMex(Fi);

    delete pm;
    break;
  }

  default:
    mexErrMsgTxt("First argument must be 'init' or 'eval'.");
  }
}
