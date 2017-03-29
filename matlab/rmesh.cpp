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
#ifdef USING_UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include "mex.h"
#include "util/matlab/MexUtil.hpp"
#include "../src/RectMeshUnstruct_pri.hpp"
#include "../src/BrickMeshBemBuilder_pri.hpp"
using namespace std;
using namespace util;
using namespace rmesh;

#include <sys/time.h>
using namespace std;

static inline void time_stamp (timeval& t) { gettimeofday(&t, 0); }

static inline double time_diff (timeval& t1, timeval& t2) {
  static const double us = 1.0e6;
  return (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
}

// Clean up command.
static bool GetCommand (const mxArray* mcmd, string& fn) {
  vector<char> vfn;
  if (!GetStringv(mcmd, vfn)) return false;
  for (int i = 0; i < vfn.size() - 1; i++) vfn[i] = tolower(vfn[i]);
  fn = string(&vfn[0]);
  return true;
}

static string _rmesh_filename = "";
#ifdef USING_UNIX
static struct stat _stats[2];
#endif
static RectMeshUnstruct* _rmu = NULL;
static RmuAnalyzer* _ra = NULL;

static bool IsSameStat (const string& fn, size_t idx) {
#ifdef USING_UNIX
  struct stat s;
  if (stat(fn.c_str(), &s)) return false;
  return s.st_mtime == _stats[idx].st_mtime;
#else
  return false;
#endif  
}

static void SetStat (const string& fn, size_t idx) {
#ifdef USING_UNIX
  stat(fn.c_str(), _stats + idx);
#endif  
}

static void FreeRa () {
  if (_ra) DeleteRmuAnalyzer(_ra);
  _ra = NULL;
}

static void FreeRmu () {
  FreeRa(); // RmuAnalyzer can't exist without RectMeshUnstruct.
  if (_rmu) DeleteRectMeshUnstruct(_rmu);
  _rmu = NULL;
  _rmesh_filename = "";
}

static void NewRa () {
  assert(_rmu);
  FreeRa();
  Boundaries b;
  _ra = NewRmuAnalyzer(_rmu, b);
}

static void LoadOrNewRa () {
  assert(_rmu);
  if (_ra && IsSameStat(_rmesh_filename + ".ra", 1)) return;
  try {
    // Use the associated .ra if it exists.
    FreeRa();
    _ra = NewRmuAnalyzer(_rmu, _rmesh_filename + ".ra");
    SetStat(_rmesh_filename + ".ra", 1);
  } catch (const FileException& e) {
    NewRa();
    // Init _stat[1] with something.
    SetStat(_rmesh_filename + ".rmu", 1);
  }
}

static void CallCubicInterp (
  const RectMeshUnstruct* rmu, const RmuAnalyzer* ra, vector<double>& z,
  const vector<double>& xi, const vector<double>& yi, vector<double>& zi)
{
  MeshAnalyzer ma(rmu, ra);
  CubicInterpolatorMatrix cim(&ma);
  zi.resize(xi.size());
#if 1
  timeval t1, t2;
  for (size_t rep = 0; rep < 2; rep++) {
    time_stamp(t1);
#endif
  for (size_t i = 0; i < xi.size(); i++) {
    TriTag tag;
    zi[i] = 0;
    if (!ma.GetTriTag(xi[i], yi[i], &tag)) continue;
    vector<double> vx(1, xi[i]), vy(1, yi[i]), vz(1);
    vector<TriTag> vtag(1, tag);
    vector<RectId> svis;
    cim.GetSupportingVertices(tag.id, svis);
    for (size_t j = 0; j < svis.size(); j++) {
      RectId id = svis[j];
      cim.Call(id, vx, vy, vtag, vz);
      zi[i] += z[id] * vz[0];
    }
  }
#if 1
    time_stamp(t2);
    printf("rep %d et %1.1f\n", rep, time_diff(t1, t2));
  }
#endif
}

static void mexExitFunction () { FreeRmu(); }

/* Implements
     id = rmesh('read', filename)
     rmesh('free')
     ids = rmesh('getids', x, y) with numel(x) == numel(y)
     rs = rmesh('getrects')
     g = rmesh('getadjgraph')
     [tri xs] = rmesh('gettri')
     tis = rmesh('gettriidxs', x, y)
     z = rmesh('linterp', z_int, z_bdy (E, N, W, S), xi, yi)
     z = rmesh('cinterp', z_int, z_bdy (E, N, W, S), xi, yi)
     z = rmesh('linterpe', z_int, xi, yi)
     z = rmesh('cinterpe', z_int, xi, yi)
     A = rmesh('cinterpm', z_bdy (E, N, W, S), xi, yi)
     A = rmesh('cinterpem', xi, yi)
 */
void mexFunction (int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
  string fn;
  if (nrhs < 1 || !GetCommand(prhs[0], fn))
    mexErrMsgTxt("[...] = rmesh('cmd',...)");
  if (fn == "read") {
    mexAtExit(&mexExitFunction);
    if (nlhs != 0 || nrhs != 2) mexErrMsgTxt("rmesh('read', filename)");
    string filename;
    if (!GetString(prhs[1], filename))
      mexErrMsgTxt("Arg 2 should be a string.");
    bool need_rmu = _rmesh_filename.empty()
      || _rmesh_filename != filename
      || !IsSameStat(_rmesh_filename + ".rmu", 0);
    if (need_rmu) {
      FreeRmu();
      try {
        _rmu = NewRectMeshUnstruct(filename + ".rmu");
        _rmesh_filename = filename;
        SetStat(_rmesh_filename + ".rmu", 0);
      } catch (const FileException& e) {
        mexErrMsgTxt(e.GetMsg().c_str());
      }
    }
  } else if (fn == "free") {
    if (nlhs != 0 || nrhs != 1) mexErrMsgTxt("rmesh('free')");
    FreeRmu();
  } else if (fn == "getids") {
    if (!_rmu) mexErrMsgTxt("No current RectMeshUnstruct.");
    if (nlhs != 1 || nrhs != 3 ||
        mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]))
      mexErrMsgTxt("ids = rmesh('getids', x, y) with numel(x) == numel(y)");
    int n = mxGetNumberOfElements(prhs[1]);
    double* px = mxGetPr(prhs[1]);
    double* py = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[1]), mxGetN(prhs[1]), mxREAL);
    double* pi = mxGetPr(plhs[0]);
    for (int i = 0; i < n; i++) {
      int id = (int) _rmu->GetRectId(px[i], py[i]);
      // Change to base-1 idx, but if it's -1, keep it.
      pi[i] = id < 0 ? id : id + 1;
    }
  } else if (fn == "getrects") {
    if (!_rmu) mexErrMsgTxt("No current RectMeshUnstruct.");
    if (nlhs != 1 || nrhs != 1) mexErrMsgTxt("rs = rmesh('getrects')");
    const vector<Rect>& rs = _rmu->GetRects();
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
    if (!_rmu) mexErrMsgTxt("No current RectMeshUnstruct.");
    if (nlhs != 2 || nrhs != 1) mexErrMsgTxt("[ag xs] = rmesh('getadjgraph')");
    if (!_ra) LoadOrNewRa();
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = _ra->GetX().Size(2);
    plhs[0] = mxCreateCellArray(2, dims);
    for (size_t i = 0; i < dims[1]; i++) {
      const vector<RectId>& nbrs = _ra->GetNbrs(i);
      mxArray* mn = mxCreateDoubleMatrix(1, nbrs.size(), mxREAL);
      double* mnp = mxGetPr(mn);
      for (size_t j = 0; j < nbrs.size(); j++) mnp[j] = (double) (nbrs[j] + 1);
      mxSetCell(plhs[0], i, mn);
    }
    plhs[1] = mxCreateDoubleMatrix(2, dims[1], mxREAL);
    memcpy(mxGetPr(plhs[1]), _ra->GetX().GetPtr(),
           2*dims[1]*sizeof(double));
  } else if (fn == "gettri") {
    if (!_rmu) mexErrMsgTxt("No current RectMeshUnstruct.");
    if (nlhs != 2 || nrhs != 1) mexErrMsgTxt("[tri xs] = rmesh('getri')");
    if (!_ra) LoadOrNewRa();
    const Matrix<RectId>& tri = _ra->GetTri();
    const RectId* ptri = tri.GetPtr();
    plhs[0] = mxCreateDoubleMatrix(3, tri.Size(2), mxREAL);
    double* pmtri = mxGetPr(plhs[0]);
    for (size_t i = 0; i < 3*tri.Size(2); i++)
      pmtri[i] = (double) (ptri[i] + 1);
    size_t nx = _ra->GetX().Size(2);
    plhs[1] = mxCreateDoubleMatrix(2, nx, mxREAL);
    memcpy(mxGetPr(plhs[1]), _ra->GetX().GetPtr(), 2*nx*sizeof(double));
  } else if (fn == "gettriidxs") {
    if (!_rmu) mexErrMsgTxt("No current RectMeshUnstruct.");
    if (nlhs != 1 || nrhs != 3 ||
        mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]))
      mexErrMsgTxt("tis = rmesh('gettriidxs', x, y)");
    if (!_ra) LoadOrNewRa();
    int n = mxGetNumberOfElements(prhs[1]);
    double* px = mxGetPr(prhs[1]);
    double* py = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[1]), mxGetN(prhs[1]), mxREAL);
    double* pi = mxGetPr(plhs[0]);
    for (int i = 0; i < n; i++) {
      // Already base-1
      TriTag tag;
      _ra->GetTriTag(px[i], py[i], &tag);
      pi[i] = (double) tag.id;
    }
  } else if (fn == "linterp" || fn == "cinterp") {
    if (!_rmu) mexErrMsgTxt("No current RectMeshUnstruct.");
    if (nlhs != 1 || nrhs != 5
        || mxGetNumberOfElements(prhs[2]) != 4
        || mxGetNumberOfElements(prhs[3]) != mxGetNumberOfElements(prhs[4])
        || mxGetNumberOfElements(prhs[1]) != _rmu->GetRects().size())
      mexErrMsgTxt(
        "z = rmesh('{l|c}interp', z_int, z_bdy (E, N, W, S), xi, yi)");
    if (!_ra) LoadOrNewRa();
    size_t ni = mxGetNumberOfElements(prhs[3]), nr = _rmu->GetRects().size();
    vector<double> xi(ni), yi(ni);
    memcpy(&xi[0], mxGetPr(prhs[3]), ni*sizeof(double));
    memcpy(&yi[0], mxGetPr(prhs[4]), ni*sizeof(double));
    vector<double> z_int(nr);
    memcpy(&z_int[0], mxGetPr(prhs[1]), nr*sizeof(double));
    vector<double> zi;
    if (fn[0] == 'l') {
      _ra->PrepInterp(z_int, mxGetPr(prhs[2]));
      _ra->LinearInterp(z_int, xi, yi, zi);
    } else {
      _ra->PrepInterp(z_int, mxGetPr(prhs[2]));
      CallCubicInterp(_rmu, _ra, z_int, xi, yi, zi);
    }
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[3]), mxGetN(prhs[3]), mxREAL);
    memcpy(mxGetPr(plhs[0]), &zi[0], ni*sizeof(double));
  } else if (fn == "linterpe" || fn == "cinterpe") {
    if (!_rmu) mexErrMsgTxt("No current RectMeshUnstruct.");
    if (nlhs != 1 || nrhs != 4
        || mxGetNumberOfElements(prhs[2]) != mxGetNumberOfElements(prhs[3])
        || mxGetNumberOfElements(prhs[1]) != _rmu->GetRects().size())
      mexErrMsgTxt("z = rmesh('{l|c}interpe', z_int, xi, yi)");
    if (!_ra) LoadOrNewRa();
    size_t ni = mxGetNumberOfElements(prhs[2]), nr = _rmu->GetRects().size();
    vector<double> xi(ni), yi(ni);
    memcpy(&xi[0], mxGetPr(prhs[2]), ni*sizeof(double));
    memcpy(&yi[0], mxGetPr(prhs[3]), ni*sizeof(double));
    vector<double> z_int(nr);
    memcpy(&z_int[0], mxGetPr(prhs[1]), nr*sizeof(double));
    vector<double> zi;
    if (fn[0] == 'l') {
      _ra->PrepExtrap(z_int);
      _ra->LinearInterp(z_int, xi, yi, zi);
    } else {
      _ra->PrepExtrap(z_int);
      CallCubicInterp(_rmu, _ra, z_int, xi, yi, zi);
    }
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[3]), mxGetN(prhs[3]), mxREAL);
    memcpy(mxGetPr(plhs[0]), &zi[0], ni*sizeof(double));
  } else {
    mexErrMsgTxt((fn + "is not a command.").c_str());
  }
}
