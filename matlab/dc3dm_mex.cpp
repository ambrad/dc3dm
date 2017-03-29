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
#include <iostream>
#include <string>
#include "mex.h"
#include "MexUtil.hpp"
#include "util/include/OpenMP.hpp"
#include "../src/RectMeshUnstruct_pri.hpp"
#include "../src/BrickMeshBemBuilder_pri.hpp"
using namespace std;
using namespace util;
using namespace rmesh;

struct Data {
  RectMeshUnstruct* rmu;
  RmuAnalyzer* ra;
  MeshAnalyzer* ma;
  CubicInterpolatorMatrix* cim;

  Data (const string& rmu_fn) throw (FileException)
    : ra(0), ma(0), cim(0), _fn(rmu_fn) {
    rmu = NewRectMeshUnstruct(_fn + ".rmu");
  }
  void InitRaIfNot () throw (FileException) {
    if (!ra) ra = NewRmuAnalyzer(rmu, _fn + ".ra");
  }
  void InitMaIfNot () throw (FileException) {
    if (!ma) {
      InitRaIfNot();
      ma = new MeshAnalyzer(rmu, ra);
    }
  }
  void InitCimIfNot () throw (FileException) {
    if (!cim) {
      InitMaIfNot();
      cim = new CubicInterpolatorMatrix(ma);
    }
  }
  ~Data () {
    if (cim) delete cim;
    if (ma) delete ma;
    if (ra) DeleteRmuAnalyzer(ra);
    DeleteRectMeshUnstruct(rmu);
  }

private:
  const string _fn;
};

vector<Data*> gds;

static void Clear (size_t id) {
  if (id < gds.size() && gds[id]) delete gds[id];
  if (id == gds.size() - 1) gds.pop_back();
  else gds[id] = 0;
}

static void ClearAll () {
  for (size_t i = 0; i < gds.size(); i++) Clear(i);
  gds.clear();
}

static size_t Init (const string& rmu_fn) {
  Data* d;
  try {
    d = new Data(rmu_fn);
  } catch (const FileException& fe) {
    mexErrMsgTxt(fe.GetMsg().c_str());
  }
  int id = -1;
  for (size_t i = 0; i < gds.size(); i++)
    if (!gds[i]) {
      id = i;
      gds[i] = d;
    }
  if (id == -1) {
    gds.push_back(d);
    id = gds.size() - 1;
  }
  return id;
}

static Data* GetData (size_t id) {
  if (id >= gds.size() || !gds[id]) mexErrMsgTxt("Invalid id.");
  return gds[id];
}

// Clean up command.
static bool GetCommand (const mxArray* mcmd, string& fn) {
  vector<char> vfn;
  if (!GetStringv(mcmd, vfn)) return false;
  for (int i = 0; i < vfn.size() - 1; i++) vfn[i] = tolower(vfn[i]);
  fn = string(&vfn[0]);
  return true;
}

static void CallCubicInterp (
  Data& d, vector<double>& z, const vector<double>& xi,
  const vector<double>& yi, vector<double>& zi)
{
  try {
    d.InitCimIfNot();
  } catch (const FileException& fe) {
    mexErrMsgTxt(fe.GetMsg().c_str());
  }
  zi.resize(xi.size());
#pragma omp parallel for
  for (size_t i = 0; i < xi.size(); i++) {
    TriTag tag;
    zi[i] = 0;
    if (!d.ma->GetTriTag(xi[i], yi[i], &tag)) continue;
    vector<double> vx(1, xi[i]), vy(1, yi[i]), vz(1);
    vector<TriTag> vtag(1, tag);
    vector<RectId> svis;
    d.cim->GetSupportingVertices(tag.id, svis);
    for (size_t j = 0; j < svis.size(); j++) {
      RectId id = svis[j];
      d.cim->Call(id, vx, vy, vtag, vz);
      zi[i] += z[id] * vz[0];
    }
  }
}

static void mexExitFunction () { ClearAll(); }

static void InitRa (Data& d) {
  try {
    d.InitRaIfNot();
  } catch (const FileException& fe) {
    mexErrMsgTxt(fe.GetMsg().c_str());
  }
}

/* Implements
     id       = dc3dm_mex('read',        filename)
                         ('free',        id)
                         ('freeall'      )
     ids      =          ('getids',      id, x, y) with numel(x) == numel(y)
     rs       =          ('getrects',    id)
     g        =          ('getadjgraph', id)
     [tri xs] =          ('gettri',      id)
     tis      =          ('gettriidxs',  id, x, y)
     z        =          ('linterp',     id, z_int, z_bdy (E, N, W, S), xi, yi)
     z        =          ('cinterp',     id, z_int, z_bdy (E, N, W, S), xi, yi)
     z        =          ('linterpe',    id, z_int, xi, yi)
     z        =          ('cinterpe',    id, z_int, xi, yi) */
void mexFunction (int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
  string fn;
  if (nrhs < 1 || !GetCommand(prhs[0], fn))
    mexErrMsgTxt("[...] = dc3dm_mex('cmd', ...)");

  if (fn == "read") {
    mexAtExit(&mexExitFunction);
    if (nlhs != 1 || nrhs != 2)
      mexErrMsgTxt("rid = dc3dm_mex('read', filename)");
    string filename;
    if (!GetString(prhs[1], filename))
      mexErrMsgTxt("Arg 2 should be a string.");
    size_t id = Init(filename);
    //printf("dc3dm_mex: read %s as %ld\n", filename.c_str(), id);
    plhs[0] = mxCreateDoubleScalar(id);
    return;
  } else if (fn == "freeall") {
    ClearAll();
    return;
  }

  if (nrhs < 2 || mxGetNumberOfElements(prhs[1]) != 1)
    mexErrMsgTxt("[...] = dc3dm_mex(command, rid, ...)");
  if (fn == "free") {
    if (nlhs != 0 || nrhs != 2) mexErrMsgTxt("dc3dm_mex('free', rid)");
    size_t rid = (size_t) mxGetScalar(prhs[1]);
    Clear(rid);
    return;
  }
  Data& d = *GetData((size_t) mxGetScalar(prhs[1]));
  nrhs -= 2;
  prhs += 2;

  if (fn == "getids") {
    if (nlhs != 1 || nrhs != 2 ||
        mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1]))
      mexErrMsgTxt("ids = dc3dm_mex('getids', rid, x, y) "
                   "with numel(x) == numel(y)");
    int n = mxGetNumberOfElements(prhs[0]);
    double* px = mxGetPr(prhs[0]);
    double* py = mxGetPr(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
    double* pi = mxGetPr(plhs[0]);
    for (int i = 0; i < n; i++) {
      int id = (int) d.rmu->GetRectId(px[i], py[i]);
      // Change to base-1 idx, but if it's -1, keep it.
      pi[i] = id < 0 ? id : id + 1;
    }
  } else if (fn == "getrects") {
    if (nlhs != 1 || nrhs != 0) mexErrMsgTxt("rs = dc3dm_mex('getrects', rid)");
    const vector<Rect>& rs = d.rmu->GetRects();
    size_t nr = rs.size();
    plhs[0] = mxCreateDoubleMatrix(4, nr, mxREAL);
    double* prs = mxGetPr(plhs[0]);
    for (size_t i = 0; i < nr; i++) {
      prs[0] = rs[i].x;
      prs[1] = rs[i].y;
      prs[2] = rs[i].dx;
      prs[3] = rs[i].dy;
      prs += 4;
    }
  } else if (fn == "getadjgraph") {
    if (nlhs != 2 || nrhs != 0)
      mexErrMsgTxt("[ag xs] = dc3dm_mex('getadjgraph', rid)");
    InitRa(d);
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = d.ra->GetX().Size(2);
    plhs[0] = mxCreateCellArray(2, dims);
    for (size_t i = 0; i < dims[1]; i++) {
      const vector<RectId>& nbrs = d.ra->GetNbrs(i);
      mxArray* mn = mxCreateDoubleMatrix(1, nbrs.size(), mxREAL);
      double* mnp = mxGetPr(mn);
      for (size_t j = 0; j < nbrs.size(); j++) mnp[j] = (double) (nbrs[j] + 1);
      mxSetCell(plhs[0], i, mn);
    }
    plhs[1] = mxCreateDoubleMatrix(2, dims[1], mxREAL);
    memcpy(mxGetPr(plhs[1]), d.ra->GetX().GetPtr(),
           2*dims[1]*sizeof(double));
  } else if (fn == "gettri") {
    if (nlhs != 2 || nrhs != 0)
      mexErrMsgTxt("[tri xs] = dc3dm_mex('getri', rid)");
    InitRa(d);
    const Matrix<RectId>& tri = d.ra->GetTri();
    const RectId* ptri = tri.GetPtr();
    plhs[0] = mxCreateDoubleMatrix(3, tri.Size(2), mxREAL);
    double* pmtri = mxGetPr(plhs[0]);
    for (size_t i = 0; i < 3*tri.Size(2); i++)
      pmtri[i] = (double) (ptri[i] + 1);
    size_t nx = d.ra->GetX().Size(2);
    plhs[1] = mxCreateDoubleMatrix(2, nx, mxREAL);
    memcpy(mxGetPr(plhs[1]), d.ra->GetX().GetPtr(), 2*nx*sizeof(double));
  } else if (fn == "gettriidxs") {
    if (nlhs != 1 || nrhs != 2 ||
        mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1]))
      mexErrMsgTxt("tis = dc3dm_mex('gettriidxs', rid, x, y)");
    InitRa(d);
    int n = mxGetNumberOfElements(prhs[0]);
    double* px = mxGetPr(prhs[0]);
    double* py = mxGetPr(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
    double* pi = mxGetPr(plhs[0]);
    for (int i = 0; i < n; i++) {
      // Already base-1
      TriTag tag;
      d.ra->GetTriTag(px[i], py[i], &tag);
      pi[i] = (double) tag.id;
    }
  } else if (fn == "linterp" || fn == "cinterp") {
    if (nlhs != 1 || nrhs != 4
        || mxGetNumberOfElements(prhs[1]) != 4
        || mxGetNumberOfElements(prhs[2]) != mxGetNumberOfElements(prhs[3])
        || mxGetNumberOfElements(prhs[0]) != d.rmu->GetRects().size())
      mexErrMsgTxt(
        "z = dc3dm_mex('{l|c}interp', rid, z_int, z_bdy (E, N, W, S), xi, yi)");
    InitRa(d);
    size_t ni = mxGetNumberOfElements(prhs[2]), nr = d.rmu->GetRects().size();
    vector<double> xi(ni), yi(ni);
    memcpy(&xi[0], mxGetPr(prhs[2]), ni*sizeof(double));
    memcpy(&yi[0], mxGetPr(prhs[3]), ni*sizeof(double));
    vector<double> z_int(nr);
    memcpy(&z_int[0], mxGetPr(prhs[0]), nr*sizeof(double));
    vector<double> zi;
    if (fn[0] == 'l') {
      d.ra->PrepInterp(z_int, mxGetPr(prhs[1]));
      d.ra->LinearInterp(z_int, xi, yi, zi);
    } else {
      d.ra->PrepInterp(z_int, mxGetPr(prhs[1]));
      CallCubicInterp(d, z_int, xi, yi, zi);
    }
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[2]), mxGetN(prhs[2]), mxREAL);
    memcpy(mxGetPr(plhs[0]), &zi[0], ni*sizeof(double));
  } else if (fn == "linterpe" || fn == "cinterpe") {
    if (nlhs != 1 || nrhs != 3
        || mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2])
        || mxGetNumberOfElements(prhs[0]) != d.rmu->GetRects().size())
      mexErrMsgTxt("z = dc3dm_mex('{l|c}interpe', rid, z_int, xi, yi)");
    InitRa(d);
    size_t ni = mxGetNumberOfElements(prhs[1]), nr = d.rmu->GetRects().size();
    vector<double> xi(ni), yi(ni);
    memcpy(&xi[0], mxGetPr(prhs[1]), ni*sizeof(double));
    memcpy(&yi[0], mxGetPr(prhs[2]), ni*sizeof(double));
    vector<double> z_int(nr);
    memcpy(&z_int[0], mxGetPr(prhs[0]), nr*sizeof(double));
    vector<double> zi;
    if (fn[0] == 'l') {
      d.ra->PrepExtrap(z_int);
      d.ra->LinearInterp(z_int, xi, yi, zi);
    } else {
      d.ra->PrepExtrap(z_int);
      CallCubicInterp(d, z_int, xi, yi, zi);
    }
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[1]), mxGetN(prhs[1]), mxREAL);
    memcpy(mxGetPr(plhs[0]), &zi[0], ni*sizeof(double));
  } else {
    mexErrMsgTxt((fn + "is not a command.").c_str());
  }
}
