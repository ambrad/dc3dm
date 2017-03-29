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
#include <map>
#include <limits>
#include <sstream>
#include "hmmvp/include/Compress.hpp"
#include "hmmvp/include/Hmat.hpp"
#include "util/include/LinAlg.hpp"
#include "util/include/Util.hpp"
#include "util/include/CodeAnalysis.hpp"
#include "util/include/IO.hpp"
#include "util/include/OpenMP.hpp"
#include "Tri2.hpp"
#include "BrickMeshBemBuilder_pri.hpp"

// I knowingly do some thread unsafe counting. But turn that off when I'm
// running this through valgrind --tool=helgrind.
#define __incr__(v) (v)++
//#define __incr__(v)
#define incr_ctr(fld) __incr__(_ctr.fld)

// If true, if a source is rsr_InterpOnly wrt the rcv, then compute the src's
// element's influence at the center of the rcv element rather than at each rcv
// subelement. [Later:] I don't like the result. The solution is only slightly
// less accurate, but the saved work is completely negligible, so I'd rather do
// it and get the slight accuracy bump.
static const bool g_rsr_InterpOnly_whole = false;

// A test idea. Not used.
//todo-remove
#ifdef SHEAR_SD_FACTOR
static const size_t g_shear_sd_factor = 4;
#endif

namespace util {
namespace rmesh {

template<typename T>
static void PrintMesh (const MeshAnalyzer& ma, const string& name,
                       const vector<T>& q) {
  const Rect& d = ma.GetDomain();
  int d_sz = ma.UnitDx(d.dx);
  printf("%s = [", name.c_str());
  for (int iy = 0; iy < d_sz; iy++) {
    double y = d.y + (iy + 0.5)*ma.uy();
    for (int ix = 0; ix < d_sz; ix++) {
      RectId id = ma.GetRectId(d.x + (ix + 0.5)*ma.ux(), y);
      cout << q[id] << " ";
    }
    if (iy == d_sz-1) printf("];\n"); else printf("\n");
  }
}

template<typename T>
static void MakeMeshArray (const MeshAnalyzer& ma, const vector<T>& q,
                           Matrix<double>& A) {
  const Rect& d = ma.GetDomain();
  int d_sz = ma.UnitDx(d.dx);
  A.Resize(d_sz, d_sz);
  for (int iy = 1; iy <= d_sz; iy++) {
    double y = d.y + (iy - 0.5)*ma.uy();
    for (int ix = 1; ix <= d_sz; ix++) {
      RectId id = ma.GetRectId(d.x + (ix - 0.5)*ma.ux(), y);
      A(iy,ix) = q[id];
    }}
}

inline ostream& operator<< (ostream& os, const Rect& r) {
  return os << "[(" << r.x << ", " << r.y << ") " << r.dx << ", " << r.dy
            << ")]";
}

RmeshBfSearcher::RmeshBfSearcher (const MeshAnalyzer* ma) : _ma(*ma) {
  // Init the mask. This is the only time an O(nx()) operation is performed. For
  // this reason, create an RmeshBfSearcher just once and call
  // FindRectIdsByDist()/FindRectsByHops() as many times as desired.  Because of
  // [D1], even with periodic BCs, nx() is an upper bound on the set size.
  _in.resize(_ma.nx(), false);
}

void RmeshBfSearcher::Clear () {
  // Clear _in.
  for (size_t i = 0; i < _rs.size(); i++) _in[_rs[i]] = false;
  for (size_t i = 0; i < _discarded.size(); i++) _in[_discarded[i]] = false;
  // Clear everything else.
  _rs.clear();
  _supports.clear();
  _discarded.clear();
}

template<typename T>
ostream& operator<< (ostream& os, const vector<T>& v) {
  for (size_t i = 0; i < v.size(); i++) os << i << " ";
  return os;
}

void RmeshBfSearcher::
FindRectIdsByDist (RectId root_id, double dist, size_t support_layers) {
  assert(support_layers <= 2);
  Clear();

  // BFS starting at root_id to get nbrs and 0 or 1 layers of supports.
  queue<RectId> to_examine;
  RectId current_id = root_id;
  _rs.push_back(root_id);
  _in[root_id] = true;
  do {
    // Examine the neighbors of current_id.
    const vector<RectId>& g_nbrs = _ma.GetNbrs(current_id);
    for (size_t i = 0; i < g_nbrs.size(); i++) {
      RectId nid = g_nbrs[i];
      if (!_in[nid]) {
        // Boundary to go through from root_id to nid, if any.
        Dir::Enum bdy_dir = _ma.GetBoundaryDir(nid);
        bool near = _ma.Distance(root_id, nid) < dist;
        if (bdy_dir == Dir::inside && near) {
          // Nearby rect.
          _rs.push_back(nid);
          to_examine.push(nid);
        } else if (near // Nearby bdy.
                   || support_layers > 0 // One hop away from a nearby id.
                   ) {
          _supports.push_back(nid);
        } else {
          _discarded.push_back(nid);
        }
        _in[nid] = true; // We've looked at it.
      }}
      
    if (to_examine.empty()) break;
    current_id = to_examine.front();
    to_examine.pop();
  } while (true);

  // Add additional layers of support.
  size_t supports_start = 0;
  for (size_t sl = 1; sl < support_layers; sl++) {
    size_t ns = _supports.size();
    for (size_t j = supports_start; j < ns; j++) {
      const vector<RectId>& g_nbrs = _ma.GetNbrs(_supports[j]);
      for (size_t i = 0; i < g_nbrs.size(); i++) {
        RectId nid = g_nbrs[i];
        if (!_in[nid]) {
          _supports.push_back(nid);
          _in[nid] = true;
        }}}

    supports_start = ns;
  }

  _rs.insert(_rs.end(), _supports.begin(), _supports.end());
}

void RmeshBfSearcher::FindRectsByHops (RectId root_id, size_t nhops) {
  Clear();
  _rs.push_back(root_id);
  if (root_id >= _ma.nr()) {
    // If root_id is a boundary point, we're going to remove it from _rs at the
    // end.
    _discarded.push_back(root_id);
  }
  _in[root_id] = true;
  size_t rs_start = 0;
  for (size_t ih = 0; ih < nhops; ih++) {
    size_t nrs = _rs.size();
    for (size_t ir = rs_start; ir < nrs; ir++) {
      const vector<RectId>& nbrs = _ma.GetNbrs(_rs[ir]);
      for (size_t in = 0; in < nbrs.size(); in++) {
        RectId nid = nbrs[in];
        if (_in[nid]) continue;
        _in[nid] = true;
        if (nid < _ma.nr()) // Want only rects.
          _rs.push_back(nid);
        else
          _discarded.push_back(nid);
      }
    }
    rs_start = nrs;
  }
  if (root_id >= _ma.nr()) _rs.erase(_rs.begin());
}

namespace {
struct NbrhdsBuilder {
  const MeshAnalyzer& ma;
  RmeshBfSearcher bfs;
  NbrhdsBuilder(const MeshAnalyzer& ima) : ma(ima), bfs(&ima) {}
};

// Set minimum neighborhood sizes.
void SetMinNhszs (NbrhdsBuilder& nb, double dist_fac, vector<size_t>& nhszs) {
  nhszs.resize(nb.ma.nr(), 0);
  if (dist_fac == 0) return;
  for (size_t i = 0; i < nhszs.size(); i++) {
    double dfncid = dist_fac*nb.ma.Nc()[i];
    size_t dfnci = dfncid > numeric_limits<size_t>::max()
      ? numeric_limits<size_t>::max() : (size_t) dfncid;
    assert (dfncid == 0 || dfnci > 0);
    nb.bfs.FindRectIdsByDist(i, dfnci, 0);
    nhszs[i] = std::max(nhszs[i], dfnci);
    const vector<RectId>& nbrs = nb.bfs.GetAll();
    for (size_t in = 0, n = nb.bfs.GetSupportsStartIdx(); in < n; in++) {
      assertpr
        (nbrs[in] < nhszs.size(),
         "nrect = %d nelem = %d nhszs.size() = %d in = %d nbrs[in] = %d",
         (int) nb.ma.nr(), (int) nb.ma.nx(),
         (int) nhszs.size(), (int) in, (int) nbrs[in]);
      nhszs[nbrs[in]] = std::max(nhszs[nbrs[in]], dfnci);
    }}
}

void BuildRcvNbrhd (NbrhdsBuilder& nb, size_t support, RcvNbrhd& rn,
                    bool dist_fac_is_0) {
  rn.nc = nb.ma.Nc()[rn.id];
  for (size_t i = 0; i < 4; i++) rn.bdy_relations[i] = rsr_Normal;
  rn.periodic_normal = true;
  if (dist_fac_is_0) return;

  nb.bfs.FindRectIdsByDist(rn.id, rn.nhsz, support);
  rn.ins = nb.bfs.GetSupport();
  std::sort(rn.ins.begin(), rn.ins.end());
  { const vector<RectId> ns = nb.bfs.GetAll();
    for (size_t i = 0; i < nb.bfs.GetSupportsStartIdx(); i++)
      rn.nc = std::min(rn.nc, nb.ma.Nc()[ns[i]]); }

#ifdef SHEAR_SD_FACTOR
  if (rn.nc == nb.ma.Nc()[rn.id]) rn.nc /= g_shear_sd_factor;
#endif

  // Now deal with v-BC boundaries. Search ins backwards to get all the boundary
  // points.
  for (int j = rn.ins.size() - 1; j >= 0; j--) {
    if (rn.ins[j] < nb.ma.nr()) break;
    Dir::Enum bdy = Dir::ToCardinal(nb.ma.GetBoundaryDir(rn.ins[j]));
    if (nb.ma.GetBoundaries().GetBC(bdy) != Boundaries::bc_velocity) continue;
    if (rn.bdy_relations[(int) bdy] != rsr_InNbrhd) {
      rn.bdy_relations[(int) bdy] =
        (nb.ma.Distance(rn.id, rn.ins[j]) < rn.nhsz)
        ? rsr_InNbrhd : rsr_InterpOnly;
    }
  }

  // If there is a layer of primary sources that are rsr_Normal w.r.t. this rcv,
  // and the layer separates this rcv from the periodic domain on the other side
  // of this layer, then all periodic src's can be treated as rsr_Normal.
  //   First, find the max distance between this element's center and the far
  // edge of a support element.
  double max_dx = 0, max_dy = 0;
  { for (size_t i = 0; i < rn.ins.size(); i++) {
      if (rn.ins[i] >= nb.ma.nr()) continue;
      double dx, dy;
      nb.ma.Distance(rn.id, rn.ins[i], &dx, &dy);
      const Rect& r = nb.ma.GetRects()[rn.ins[i]];
      max_dx = max(max_dx, dx + r.dx);
      max_dy = max(max_dy, dy + r.dy);
    }
    const Rect& r = nb.ma.GetRects()[rn.id];
    max_dx += r.dx;
    max_dy += r.dy; }
  rn.periodic_normal =
    (nb.ma.GetBoundaries().GetBC(Dir::E) != Boundaries::bc_periodic
     || nb.ma.UnitDx(2*max_dx) < nb.ma.UnitDx(nb.ma.GetDomain().dx)) &&
    (nb.ma.GetBoundaries().GetBC(Dir::N) != Boundaries::bc_periodic
     || nb.ma.UnitDy(2*max_dy) < nb.ma.UnitDy(nb.ma.GetDomain().dy));
#if 0
  if (!rn.periodic_normal)
    printf("rcv %ld 2*max_dx %d 2*max_dy %d domain.dx %d domain.dy %d\n",
           rn.id, nb.ma.UnitDx(2*max_dx), nb.ma.UnitDy(2*max_dy),
           nb.ma.UnitDx(nb.ma.GetDomain().dx),
           nb.ma.UnitDy(nb.ma.GetDomain().dy));
#endif
}

void BuildSrcNbrhd (NbrhdsBuilder& nb, size_t support, SrcNbrhd& sn,
                    bool dist_fac_is_0) {
  if (!dist_fac_is_0) {
    nb.bfs.FindRectsByHops(sn.id, support);
    sn.ns = nb.bfs.GetAll();
  }
}

void BuildNbrhds (const MeshAnalyzer& ma, size_t support, double dist_fac,
                  vector<RcvNbrhd>& rns, vector<SrcNbrhd>& sns) {
  NbrhdsBuilder nb(ma);
  const size_t nr = ma.nr();
  const size_t nx = ma.nx();
  // Work out neighborhood sizes.
  vector<size_t> nhszs;
  SetMinNhszs(nb, dist_fac, nhszs);
  // Receiver rects.
  rns.resize(nr);
  for (size_t i = 0; i < rns.size(); i++) {
    rns[i].id = (RectId) i;
    rns[i].nhsz = nhszs[i];
    BuildRcvNbrhd(nb, support, rns[i], dist_fac == 0);
  }
  // Source rects and boundary points.
  sns.resize(nx);
  for (size_t i = 0; i < sns.size(); i++) {
    sns[i].id = (RectId) i;
    BuildSrcNbrhd(nb, support, sns[i], dist_fac == 0);
  }
}
}

RcvSrcRelation BrickMeshBemBuilder::
GetRelation (const RcvNbrhd& rn, RectId src_id) const {
  if (src_id < _nx) {
    if (src_id < _nr
        && (rn.id == src_id || _ma.Distance(rn.id, src_id) < rn.nhsz))
      return rsr_InNbrhd;
    if (std::binary_search(rn.ins.begin(), rn.ins.end(), src_id))
      return rsr_InterpOnly;
    return rsr_Normal;
  } else {
    // We're dealing with a v-BC boundary rect.
    return rn.bdy_relations[src_id - _nx];
  }
}

BrickMeshBemBuilder::
BrickMeshBemBuilder (
  const GreensFn* gf, const MeshAnalyzer* ma, const InterpolatorMatrix* im,
  double dist_fac, size_t n_per_layers)
  : _gf(*gf), _ma(*ma), _im(*im), _n_per_layers(n_per_layers),
    _full_quality(true), _use_lra(true), _dist_fac_is_0(dist_fac == 0)
{
  _nr = _ma.nr();
  _nx = _ma.nx();
  SetOmpNthreads(1);

  int np = (int) (_ma.GetBoundaries().GetBC(Dir::E) == Boundaries::bc_periodic)
    + (int) (_ma.GetBoundaries().GetBC(Dir::N) == Boundaries::bc_periodic);
  _any_periodic = np > 0;

#ifdef SHEAR_SD_FACTOR
  _ma.ResetUnitFactor(g_shear_sd_factor);
#endif
  BuildNbrhds(_ma, _im.GetNbrTriLayersInSupport(), dist_fac, _rns, _sns);
  BuildBdyData();
}

BrickMeshBemBuilder::
BrickMeshBemBuilder (const GreensFn* gf, const MeshAnalyzer* ma,
                     const InterpolatorMatrix* im, const string& filename)
  throw (FileException)
  : _gf(*gf), _ma(*ma), _im(*im), _full_quality(true), _use_lra(true)
{
  FILE* fid = fopen(filename.c_str(), "r");
  if (!fid) throw FileException("Can't read " + filename);
  Deserialize(fid);
  fclose(fid);

  if (_nr != _ma.nr() || _nx != _ma.nx())
    throw FileException("_im.nr|x() != _nr|x");

  SetOmpNthreads(1);
  
  _dist_fac_is_0 = true;
  for (size_t i = 0; i < _rns.size(); i++)
    if (_rns[i].nhsz > 0) {
      _dist_fac_is_0 = false;
      break;
    }

#ifdef SHEAR_SD_FACTOR
  _ma.ResetUnitFactor(g_shear_sd_factor);
#endif
}

BrickMeshBemBuilder::~BrickMeshBemBuilder () {}

void BrickMeshBemBuilder::Serialize (const string& filename) const
  throw (FileException) {
  FILE* fid = fopen(filename.c_str(), "w");
  if (!fid) throw FileException("Can't read " + filename);
  Serialize(fid);
  fclose(fid);
}

namespace ser {
  void Serialize (FILE* fid, const RcvNbrhd* rn) {
    write(&rn->id, 1, fid);
    write(&rn->nhsz, 1, fid);
    Write(rn->ins, fid);
    write(rn->bdy_relations, 4, fid);
    write(&rn->nc, 1, fid);
    write(&rn->periodic_normal, 1, fid);
  }

  void Deserialize (FILE* fid, RcvNbrhd* rn) {
    read(&rn->id, 1, fid);
    read(&rn->nhsz, 1, fid);
    Read(rn->ins, fid);
    read(rn->bdy_relations, 4, fid);
    read(&rn->nc, 1, fid);
    read(&rn->periodic_normal, 1, fid);
  }

  void Serialize (FILE* fid, const SrcNbrhd* sn) {
    write(&sn->id, 1, fid);
    Write(sn->ns, fid);
    Write(sn->fbps, fid);
    Write(sn->fbps_wts, fid);
  }

  void Deserialize (FILE* fid, SrcNbrhd* sn) {
    read(&sn->id, 1, fid);
    Read(sn->ns, fid);
    Read(sn->fbps, fid);
    Read(sn->fbps_wts, fid);
  }

  void Serialize (FILE* fid, const BdyData* bd) {
    write(&bd->id, 1, fid);
    Write(bd->bdy_pts, fid);
    write(&bd->r.x, 4, fid);
  }

  void Deserialize (FILE* fid, BdyData* bd) {
    read(&bd->id, 1, fid);
    Read(bd->bdy_pts, fid);
    read(&bd->r.x, 4, fid);
  }
}

void BrickMeshBemBuilder::Serialize (FILE* fid) const throw (FileException) {
  size_t n;

  write("bmb ", 4, fid);
  write(&_any_periodic, 1, fid);
  write(&_nr, 1, fid);
  write(&_nx, 1, fid);
  
  n = _rns.size(); write(&n, 1, fid);
  for (size_t i = 0; i < n; i++) ser::Serialize(fid, &_rns[i]);
  
  n = _sns.size(); write(&n, 1, fid);
  for (size_t i = 0; i < n; i++) ser::Serialize(fid, &_sns[i]);

  for(size_t i = 0; i < 4; i++) ser::Serialize(fid, _bdys+i);
  write(&_n_per_layers, 1, fid);
}

void BrickMeshBemBuilder::Deserialize (FILE* fid) {
  size_t n;

  char tag[4];
  read(tag, 4, fid);
  if (string(tag, 4) != string("bmb ")) {
    printf("tag = %s\n", string(tag, 4).c_str());
    throw FileException("Not a BrickMeshBemBuilder serial file.");
  }

  read(&_any_periodic, 1, fid);
  read(&_nr, 1, fid);
  read(&_nx, 1, fid);

  read(&n, 1, fid);
  _rns.resize(n);
  for (size_t i = 0; i < n; i++) ser::Deserialize(fid, &_rns[i]);
  
  read(&n, 1, fid);
  _sns.resize(n);
  for (size_t i = 0; i < n; i++) ser::Deserialize(fid, &_sns[i]);

  for(size_t i = 0; i < 4; i++) ser::Deserialize(fid, _bdys+i);
  read(&_n_per_layers, 1, fid);
}

size_t BrickMeshBemBuilder::GetCtrNgfs () const { return _ctr.ngfs; }
size_t BrickMeshBemBuilder::GetCtrNrgfs () const { return _ctr.nrgfs; }
size_t BrickMeshBemBuilder::GetCtrNcalls () const { return _ctr.ncalls; }

template<typename T> static inline T Square (const T& v) { return v * v; }

void BrickMeshBemBuilder::BuildBdyData () {
  // Fill in data for v-BC boundaries.
  for (size_t i = 0; i < 4; i++) {
    _bdys[i].id = _nx + i;
    Dir::Enum dir = (Dir::Enum) i;
    if (_ma.GetBoundaries().GetBC(dir) == Boundaries::bc_velocity)
      _bdys[i].r = _ma.GetBoundaries().GetRect(dir);
  }
  if (_dist_fac_is_0) return;
  // Collect bdy points.
  for (size_t i = _nr; i < _nx; i++) {
    Dir::Enum dir = _ma.GetBoundaryDir(i);
    Boundaries::BC bc;
    assert (dir != Dir::inside);
    if (Dir::IsCardinal(dir)) {
      _bdys[(int) dir].bdy_pts.push_back(i);
      bc = _ma.GetBoundaries().GetBC(dir);
    } else { // Corner point.
      Dir::Enum dir1, dir2;
      Dir::BreakUp(dir, dir1, dir2);
      const Boundaries::BC bc1 = _ma.GetBoundaries().GetBC(dir1);
      const Boundaries::BC bc2 = _ma.GetBoundaries().GetBC(dir2);
      // [c1] The following if list sets priorities on boundary conditions from
      // velocity (considered first) to free (last). This ordering matters at
      // corners.
      size_t idx;
#define do1 { bc = bc1; idx = (int) dir1; }
#define do2 { bc = bc2; idx = (int) dir2; }
      if      (bc1 == Boundaries::bc_velocity)  do1
      else if (bc2 == Boundaries::bc_velocity)  do2
      else if (bc1 == Boundaries::bc_0velocity) do1
      else if (bc2 == Boundaries::bc_0velocity) do2
      else if (bc1 == Boundaries::bc_free)      do1
      else if (bc2 == Boundaries::bc_free)      do2
      else assert (false);
      _bdys[idx].bdy_pts.push_back(i);
#undef do1
#undef do2

#ifndef NDEBUG
      { Dir::Enum side;
        Boundaries::BC bc_;
        _ma.GetCornerBC(dir, side, bc_);
        assert ((int) side < 4);
        assert (side == (Dir::Enum) idx);
        assert (bc_ == bc); }
#endif
    }
  }
  // If a bdy point belongs to a free surface BC, we have to say how to
  // extrapolate the value at that point. The fbps(_bdy) fields contain the data
  // to do that.
  for (size_t k = 0; k < 4; k++) {
    if (_ma.GetBoundaries().GetBC((Dir::Enum) k) != Boundaries::bc_free)
      continue;
    BdyData& bd = _bdys[k];
    for (size_t i = 0; i < bd.bdy_pts.size(); i++) {
      RectId id = bd.bdy_pts[i];
      const SrcNbrhd& sn = _sns[id];
      // True because it's obviously true on the interior boundary (since one
      // bdy vertex is shared between two rects), and [c1] makes it true on the
      // corners (since a v-free corner point is v, and a periodic-free corner
      // point gets the two periodically adjacent rects).
      //   If we're using cubic interp, then the first two are the adjacent
      // rects, and these are the only two we want.
      assert ((_im.GetNbrTriLayersInSupport() == 1 && sn.ns.size() == 2) ||
              (_im.GetNbrTriLayersInSupport() == 2 && sn.ns.size() >= 2));
      assert (sn.ns[0] < _nr && sn.ns[1] < _nr);
      if (sn.ns.size() > 2)
        CalcFbpsWts(sn);
      else {
        double areas[2];
        for (size_t j = 0; j < 2; j++) {
          _sns[sn.ns[j]].fbps.push_back(id);
          const Rect& r = _ma.GetRects()[sn.ns[j]];
          areas[j] = r.dx * r.dy;
        }
        for (size_t j = 0; j < 2; j++)
          _sns[sn.ns[j]].fbps_wts.push_back(areas[j] / (areas[0] + areas[1]));
      }
    }
  }
}

//todo-opt I don't actually need Q in the QR factorization.
/* Fit a plane to (xi, yi) points and then eval at (x, y).
       z = c + m1 x + m2 y
       A = [ones xi(:) yi(:)]
       z = [1 x y] inv(A' A) A' b = w' b.
   b is the vector of (unknown) function values at each (xi, yi). Then the
   weights on each point are
       w' = [1 x y] inv(A' A) A',
   obtained by
       A = Q R;
       solve (R' R) alpha = [1 x y]';
       set w' = alpha' A'. */
void BrickMeshBemBuilder::CalcFbpsWts (const SrcNbrhd& sn) {
  const size_t n = sn.ns.size();

  Matrix<double> A(n, 3), R(3, 3);
  { const double* bx = &_ma.GetX()(1, sn.id+1);
    Matrix<double> Q(n, 3);
    for (size_t j = 0; j < n; j++) {
      assert(sn.ns[j] < _nr);
      const Rect sr = GetPrimarySrcRect(sn.id, sn.ns[j]);
      A(j+1, 1) = 1;
      A(j+1, 2) = sr.x - bx[0] + 0.5*sr.dx;
      A(j+1, 3) = sr.y - bx[1] + 0.5*sr.dy;
    }
    SkinnyQr(A, Q, R); }
  
  double alpha[3];
  { alpha[0] = 1; alpha[1] = 0; alpha[2] = 0;
    const blas_int n = 3, nrhs = 1;
    blas_int info;
    trtrs('u', 't', 'n', n, nrhs, R.GetPtr(), n, alpha, n, info);
    assert (info == 0);
    trtrs('u', 'n', 'n', n, nrhs, R.GetPtr(), n, alpha, n, info);
    assert (info == 0); }

  const double* pA = A.GetPtr();
  for (size_t j = 0; j < n; j++) {
    _sns[sn.ns[j]].fbps.push_back(sn.id);
    _sns[sn.ns[j]].fbps_wts.push_back(
        alpha[0]*pA[j] + alpha[1]*pA[n + j] + alpha[2]*pA[2*n + j]);
  }

#if 0
  // Compare with RmuAnalyzer::Extrap debug output, the routine I tested the
  // method with.
  Matrix<double> w(n);
  for (size_t j = 0; j < n; j++)
    w(j+1) = alpha[0]*pA[j] + alpha[1]*pA[n + j] + alpha[2]*pA[2*n + j];
  cout << "bws{" << sn.id+1 << "} = " << w << ";" << endl;
#endif
}

// Get the possibly periodic source Rect nearest the receiver.
inline Rect BrickMeshBemBuilder::
GetPrimarySrcRect (RectId rcv, RectId src, Dir::Enum* idir) const {
  if (src < _nr) {
    if (!_any_periodic) return _ma.GetRects()[src];
    Dir::Enum dir;
    _ma.Distance(rcv, src, &dir);
    Rect r = _ma.GetRects()[src];
    _ma.TranslateRect(r, dir);
    if (idir) *idir = dir;
    return r;
  } else {
    assert (src >= _nx && src < _nx + 4);
    return _bdys[src - _nx].r;
  }
}

namespace bdy {
// The implementation of a v BC implicitly is adding a bunch of columns of the
// interpolation matrix together. Each such column corresponds to a boundary
// point. The Accumulator manages this operation.

struct Key {
  double xi, yi;
  Rect r;
  TriTag in;

  Key(double ixi, double iyi, const Rect& ir, const TriTag& ini)
    : xi(ixi), yi(iyi), r(ir), in(ini) {}
};

inline bool operator< (const Key& k1, const Key& k2) {
  if (k1.xi   < k2.xi  ) return true;
  if (k1.xi   > k2.xi  ) return false;
  if (k1.yi   < k2.yi  ) return true;
  if (k1.yi   > k2.yi  ) return false;
  // I realized I don't need this next bit. Source centers are enough to
  // differentiate them.
  assert (k1.r.dx == k2.r.dx && k1.r.dy == k2.r.dy);
#if 0
  if (k1.r.dx < k2.r.dx) return true;
  if (k1.r.dx > k2.r.dx) return false;
  if (k1.r.dy < k2.r.dy) return true;
#endif
  return false;
}

class Accumulator {
public:
  void Accumulate(const BrickMeshBemBuilder::Work& w);
  void Finalize(BrickMeshBemBuilder::Work& w) const;
private:
  map<Key, double> _ks;
};

void Accumulator::Accumulate (const BrickMeshBemBuilder::Work& w) {
  for (size_t i = 0; i < w.xi.size(); i++) {
    pair<Key, double> p(Key(w.xi[i], w.yi[i], w.srcs[i], w.in[i]), w.zi[i]);
    map<Key, double>::iterator it = _ks.find(p.first);
    if (it == _ks.end()) _ks.insert(p);
    else it->second += w.zi[i];
  }
}

void Accumulator::Finalize (BrickMeshBemBuilder::Work& w) const {
  // Load work with what we found. xi, yi, in aren't actually used
  // subsequently, so maybe I'll drop them from here at some point.
  w.ClearSrcs();
  for (map<Key, double>::const_iterator it = _ks.begin(); it != _ks.end(); ++it)
    if (it->second != 0) {
      w.xi.push_back(it->first.xi); w.yi.push_back(it->first.yi);
      w.in.push_back(it->first.in);
      w.srcs.push_back(it->first.r);
      w.zi.push_back(it->second);
    }
}
}

inline void BrickMeshBemBuilder::Work::ClearSrcs () {
  xi.clear(); yi.clear(); zi.clear(); srcs.clear(); in.clear();
}

bool BrickMeshBemBuilder::
GetWeightedSourcesSni (const RcvNbrhd& rn, const SrcNbrhd& sn,
                       const SrcNbrhd& sni, Work& w) const {
  bool need_interp = false;
  RcvSrcRelation rsr = GetRelation(rn, sni.id);
  switch (rsr) {
  case rsr_Normal:
  case rsr_InterpOnly: {
    if (!g_rsr_InterpOnly_whole && sni.id == sn.id) {
      assert (rsr == rsr_InterpOnly);
      // I'm the source of interest, so I'm definitely a non-0 entry (and, in
      // fact, my zi is 1).
      const Rect sr = GetPrimarySrcRect(rn.id, sn.id);
      const double x = sr.x + 0.5 * sr.dx;
      const double y = sr.y + 0.5 * sr.dy;
      w.srcs.push_back(sr);
      w.xi.push_back(x); w.yi.push_back(y);
      w.in.push_back(TriTag()); _ma.GetTriTag(x, y, &w.in.back());
      w.zi.push_back(1);
    }
    // But otherwise my zi value is 0, so skip me.
  } break;
  case rsr_InNbrhd: {
    // Break the source into subrects.
    const size_t nsub = _ma.Nc()[sni.id] / rn.nc;
    assert (nsub * rn.nc == _ma.Nc()[sni.id]);
    const Rect r = GetPrimarySrcRect(rn.id, sni.id);
    const double dx = r.dx / nsub, dxh = 0.5 * dx;
    const double dy = r.dy / nsub, dyh = 0.5 * dy;
    Rect sr(0.0, r.y, dx, dy);
    for (size_t iy = 0; iy < nsub; iy++) {
      const double y = sr.y + dyh;
      sr.x = r.x;
      for (size_t ix = 0; ix < nsub; ix++) {
        const double x = sr.x + dxh;
        w.srcs.push_back(sr);
        w.xi.push_back(x); w.yi.push_back(y);
        w.in.push_back(TriTag()); _ma.GetTriTag(x, y, &w.in.back());
        w.zi.push_back(0.0);
        sr.x += dx; }
      sr.y += dy;
    }
    if (nsub > 1) need_interp = true;
    else if (sni.id == sn.id) {
      // I'm the source of interest, and I wasn't broken up, so my zi is 1. Set
      // it in case !need_interp.
      w.zi.back() = 1;
    }
  } break;
  }
  return need_interp;
}

void BrickMeshBemBuilder::
GetWeightedSources (const RcvNbrhd& rn, const SrcNbrhd& sn) const {
  const int tid = omp_get_thread_num();
  _w[tid].ClearSrcs();
  // If no source rect is broken up, then interpolation will give a single 1
  // (for sn.id) and the rest 0. So detect whether we need to call
  // _im.Call. An Accumulator is not necessary because all points are distinct
  // by construction.
  bool need_interp = false;
  for (size_t is = 0; is < sn.ns.size(); is++) {
    const SrcNbrhd& sni = _sns[sn.ns[is]];
    assert (sni.id < _nr);
    need_interp = GetWeightedSourcesSni(rn, sn, sni, _w[tid]) || need_interp;
  }
  if (need_interp)
    _im.Call(sn.id, _w[tid].xi, _w[tid].yi, _w[tid].in, _w[tid].zi);
  // Handle the free BC if it exists.
  if (!sn.fbps.empty()) {
    bdy::Accumulator ba;
    ba.Accumulate(_w[tid]);
    for (size_t fi = 0; fi < sn.fbps.size(); fi++) {
      // Collect all sources just as before. We can't use need_interp in this
      // case because the value-1 node is fid, not sn.id.
      _w[tid].ClearSrcs();
      const RectId fid = sn.fbps[fi];
      const SrcNbrhd& snf = _sns[fid];
      assert (sn.id < _nr && snf.id >= _nr);
      for (size_t is = 0; is < snf.ns.size(); is++) {
        const SrcNbrhd& sni = _sns[snf.ns[is]];
        GetWeightedSourcesSni(rn, sn, sni, _w[tid]);
      }
      _im.Call(fid, _w[tid].xi, _w[tid].yi, _w[tid].in, _w[tid].zi);
      // Weight these new sources by sn's contribution to the boundary point.
      for (size_t i = 0; i < _w[tid].zi.size(); i++)
        _w[tid].zi[i] *= sn.fbps_wts[fi];
      ba.Accumulate(_w[tid]);
    }
    ba.Finalize(_w[tid]);
  }
}

void BrickMeshBemBuilder::
GetWeightedSources (const RcvNbrhd& rn, RectId src) const {
  if (src < _nr) {
    GetWeightedSources(rn, _sns[src]);
    return;
  }
  // We're dealing with a v-BC rect.
  assert (src >= _nx);
  const BdyData& bd = _bdys[src - _nx];
  bdy::Accumulator ba;
  const int tid = omp_get_thread_num();
  for (size_t i = 0; i < bd.bdy_pts.size(); i++) {
    GetWeightedSources(rn, _sns[bd.bdy_pts[i]]);
    ba.Accumulate(_w[tid]);
  }
  ba.Finalize(_w[tid]);
  // Append the v-BC rect itself.
  _w[tid].xi.push_back(0); _w[tid].yi.push_back(0); // Not used.
  _w[tid].in.push_back(TriTag()); // Not used.
  if (g_rsr_InterpOnly_whole && GetRelation(rn, src) != rsr_InNbrhd)
    _w[tid].zi.push_back(0);
  else
    _w[tid].zi.push_back(1); // It contributes its full weight.
  _w[tid].srcs.push_back(bd.r);
}

inline double BrickMeshBemBuilder::
CallGreensFn (const Rect& sr, double rx, double ry, bool use_memo) const {
  double B;
  if (use_memo) {
    const int tid = omp_get_thread_num();
    pair<Memoizer::iterator, bool> ret = _w[tid].memo.Lookup(sr, rx, ry, &B);
    if (ret.second) {
      incr_ctr(ngfs);
      B = _gf.Call(sr, rx, ry);
      _w[tid].memo.Insert(ret, B);
    }
  } else {
    incr_ctr(ngfs);
    B = _gf.Call(sr, rx, ry);
  }
  return B;
}

inline double BrickMeshBemBuilder::
Call_const (RectId rcv, RectId src, double tol) const {
  if (src >= _nr) src += (_nx - _nr);

  incr_ctr(ncalls);
  double B = 0;
  // For generality, a same-same calculation is rsr_InNbrhd even if
  // _dist_fac_is_0, but here I want to break that rule so that if
  // _dist_fac_is_0, we don't have to build any neighborhood information. In
  // particular, SrcNbrhd::ns.empty() even though it should have its own source.
  RcvSrcRelation rsr = _dist_fac_is_0 ? rsr_Normal :
    GetRelation(_rns[rcv], src);
  if (!_full_quality || rsr == rsr_Normal ||
      (g_rsr_InterpOnly_whole && rsr == rsr_InterpOnly)) {
    const Rect psr = GetPrimarySrcRect(rcv, src);
    B = CallGreensFn(psr, _ma.GetX()(1, rcv+1), _ma.GetX()(2, rcv+1), false);
    // If the primary source rect is rsr_Normal, then certainly so are its
    // periodic images.
    if (_any_periodic
        // So we don't translate the v-BC slab.
        && src < _nr)
      B += CallOnPeriodicImages(psr, _ma.GetX()(1, rcv+1),
                                _ma.GetX()(2, rcv+1), false);
  }
  if (_full_quality && rsr != rsr_Normal) {
    const RcvNbrhd rn = _rns[rcv];
    const size_t nsub = _ma.Nc()[rcv] / rn.nc;
    assert (nsub * rn.nc == _ma.Nc()[rcv]);
    const int tid = omp_get_thread_num();
    assert (tid >= 0 && tid < (int) _w.size());
    GetWeightedSources(rn, src);

    // If we are permitted and it makes sense, use an H-matrix approximation to
    // this subelement matrix rather than computing it in full.
    static const size_t nz_thresh = 8;
    size_t nzi = 0;
    for (size_t is = 0; is < _w[tid].zi.size(); is++) {
      if (_w[tid].zi[is] != 0) nzi++;
      if (nzi == nz_thresh) break;
    }

    if (_use_lra && tol > 0 && nsub >= nz_thresh && nzi >= nz_thresh)
      B += CallHmmvp(rcv, src, tol);
    else {
      const Rect& r = _ma.GetRects()[rcv];
      const double dx = r.dx / nsub, dy = r.dy / nsub;
      double B1 = 0;
      for (size_t is = 0; is < _w[tid].zi.size(); is++) { // For each source
        if (_w[tid].zi[is] == 0) continue;
        double ry = r.y + 0.5*dy;
        for (size_t iy = 0; iy < nsub; iy++) { // For each rcv subrect
          double rx = r.x + 0.5*dx;
          for (size_t ix = 0; ix < nsub; ix++) {
            B1 += _w[tid].zi[is] * CallGreensFn(_w[tid].srcs[is], rx, ry);
            if (_any_periodic && !rn.periodic_normal
                // If src >= _nr, the final zi is for the v-BC slab, and this
                // should not be periodically repeated.
                && (src < _nr || is != _w[tid].zi.size() - 1)) {
              B1 += _w[tid].zi[is] *
                CallOnPeriodicImages(_w[tid].srcs[is], rx, ry);
            }
            rx += dx; }
          ry += dy; }
      }
      B1 /= (nsub * nsub);
      B += B1;
    }

    // Periodic images are rsr_Normal, so handle that here.
    if (_any_periodic && rn.periodic_normal && src < _nr)
      B += CallOnPeriodicImages(
        GetPrimarySrcRect(rcv, src), _ma.GetX()(1, rcv+1), _ma.GetX()(2, rcv+1),
        false);
  }

  return B;
}

class HmmvpBlockGreensFn : public hmmvp::GreensFn {
public:
  vector<double> rxs, rys;
  int _tid_check;

  HmmvpBlockGreensFn (
    const BrickMeshBemBuilder* bmb, const vector<size_t>& wis, RectId src,
    RectId rcv, const Rect* r_rcv, size_t rnsub, int tid_check)
    : _bmb(bmb), _wis(wis), _src(src), _rcv(rcv), _r(r_rcv), _rnsub(rnsub),
      _dx(r_rcv->dx / rnsub), _dy(r_rcv->dy / rnsub), _tid_check(tid_check)
  {
    const size_t nsub2 = rnsub*rnsub;
    rxs.resize(nsub2);
    rys.resize(nsub2);
    double rx = _r->x + 0.5*_dx;
    for (size_t ix = 0, k = 0; ix < rnsub; ix++) {
      double ry = _r->y + 0.5*_dy;
      for (size_t iy = 0; iy < rnsub; iy++, k++) {
        rxs[k] = rx;
        rys[k] = ry;
        ry += _dy;
      }
      rx += _dx;
    }
  }

  virtual bool
  Call (const hmmvp::CompressBlockInfo& cbi, const vector<UInt>& rs,
        const vector<UInt>& cs, double* B) const {
    //todo Need to decide if I should use memoization here.
    const bool use_memo = true;
    const int tid = omp_get_thread_num();
    // Make sure we're handling the various OpenMP parallelizations correctly.
    assert(tid == _tid_check);
    for (size_t ic = 0, k = 0; ic < cs.size(); ic++) {
      size_t is = _wis[cs[ic]-1];
      assert (_bmb->_w[tid].zi[is] != 0);
      for (size_t ir = 0; ir < rs.size(); ir++, k++) {
        const double rx = rxs[rs[ir]-1], ry = rys[rs[ir]-1];
        B[k] = _bmb->CallGreensFn(_bmb->_w[tid].srcs[is], rx, ry, use_memo);
        if (_bmb->_any_periodic && !_bmb->_rns[_rcv].periodic_normal
            && (_src < _bmb->_nr || cs[ic] != _bmb->_w[tid].zi.size()))
          B[k] += _bmb->CallOnPeriodicImages(_bmb->_w[tid].srcs[is], rx, ry,
                                             use_memo);
      }
    }
    return true;
  }

private:
  const BrickMeshBemBuilder* _bmb;
  const vector<size_t>& _wis;
  RectId _src, _rcv;
  const Rect* _r;
  size_t _rnsub;
  double _dx, _dy;
};

inline double BrickMeshBemBuilder::
CallHmmvp (RectId rcv, RectId src, double tolp) const {
  const int tid = omp_get_thread_num();
  assert (tid >= 0 && tid < (int) _w.size());
  const RcvNbrhd& rn = _rns[rcv];
  // Already called:
  //   GetWeightedSources(rn, src);
  const Rect& r = _ma.GetRects()[rcv];
  const size_t nsub = _ma.Nc()[rcv] / rn.nc, nsub2 = nsub*nsub;

  vector<size_t> wis; // _w[tid].zi[wis[i]] != 0
  wis.reserve(_w[tid].zi.size());
  for (size_t is = 0; is < _w[tid].zi.size(); is++)
    if (_w[tid].zi[is] != 0) wis.push_back(is);

  hmmvp::Hmat* hm;
  { Matrix<double> D(3, wis.size()), R(3, nsub2);
    for (size_t i = 0; i < wis.size(); i++) {
      const Rect& r = _w[tid].srcs[wis[i]];
      D(1, i+1) = r.x + 0.5*r.dx;
      D(2, i+1) = r.y + 0.5*r.dy;
      D(3, i+1) = 0;
    }
    HmmvpBlockGreensFn hmgf(this, wis, src, rcv, &r, nsub, tid);
    for (size_t i = 0; i < nsub2; i++) {
      R(1, i+1) = hmgf.rxs[i];
      R(2, i+1) = hmgf.rys[i];
      R(3, i+1) = 0;
    }
    // We don't need to worry about a periodic domain here because the source
    // elements are translated periodically to be nearest the receiver elements.
    hmmvp::Hd* hd = hmmvp::NewHd(D, R);
    hmmvp::Compressor* c = hmmvp::NewCompressor(hd, &hmgf);
    /* We want
           x' E y <= abs(x)' E abs(y)' = sum(abs(x)) sum(abs(y)) E11 <= tolp
           E11 <= tolp / (sum(abs(x)) sum(abs(y))).
       c->SetTol() with tm_mrem_abs is specified as a tolerance on
          norm(E, fro) = sqrt(m n E11^2) = sqrt(m n) E11 <= tol
          E11 <= tol / sqrt(m n),
       and so
          tol = tolp sqrt(m n) / (sum(abs(x)) sum(abs(y))).
       now, x = e / nsub^2, and so sum(abs(x)) = 1; and y = w, the vector of
       weights. */
    { double sum_w = 0;
      for (size_t i = 0; i < wis.size(); i++) sum_w += fabs(_w[tid].zi[wis[i]]);
      double tol = tolp * nsub * sqrt(wis.size()) / sum_w;
      c->SetTolMethod(hmmvp::Compressor::tm_mrem_abs);
      c->SetTol(tol); }
    c->SetOutputLevel(0);
    c->AvoidRedundantGfCalls(true);
    c->UseCompressQr(false);
    hm = c->CompressInMemory(1, 1);
    DeleteCompressor(c);
    DeleteHd(hd); }

  // 1/nsub^2 e' H w
  double b;
  { Matrix<double> w(wis.size()), Hw(nsub2);
    double* pw = w.GetPtr();
    for (size_t i = 0; i < wis.size(); i++) pw[i] = _w[tid].zi[wis[i]];
    hm->Mvp(pw, Hw.GetPtr(), 1);
    const double* pHw = Hw.GetPtr();
    b = 0;
    for (size_t i = 0, n = nsub2; i < n; i++) b += pHw[i];
    b /= (nsub2); }

  hmmvp::DeleteHmat(hm); 

  return b;
}

double BrickMeshBemBuilder::
Call (RectId rcv, RectId src, double tol) { return Call_const(rcv, src, tol); }

namespace perimg {
  class Caller {
  public:
    Caller(const BrickMeshBemBuilder* bmb, const MeshAnalyzer& ma,
           double rx, double ry, const Rect& psr, bool use_memo)
      : _bmb(bmb), _ma(ma), _psr(psr), _B(0), _rx(rx), _ry(ry),
        _use_memo(use_memo) {}

    double B () const { return _B; }

    void Call (int Dx, int Dy) {
      Rect r = _psr;
      _ma.TranslateRect(r, Dx, Dy);
      _B += _bmb->CallGreensFn(r, _rx, _ry, _use_memo);
    }
    
  private:
    const BrickMeshBemBuilder* _bmb;
    const MeshAnalyzer& _ma;
    Rect _psr;
    double _B;
    double _rx, _ry;
    bool _use_memo;
  };
}

double BrickMeshBemBuilder::
CallOnPeriodicImages (const Rect& psr, double rx, double ry,
                      bool use_memo) const {
  assert (_any_periodic);
  if (_n_per_layers == 0) return 0;
  perimg::Caller pc(this, _ma, rx, ry, psr, use_memo);
  const bool pE = _ma.GetBoundaries().GetBC(Dir::E) == Boundaries::bc_periodic;
  const bool pN = _ma.GetBoundaries().GetBC(Dir::N) == Boundaries::bc_periodic;
  for (int layer = 1; layer <= (int) _n_per_layers; layer++) {
    if (pE && pN) {
      for (int Dx = -layer; Dx <= layer; Dx++) {
        pc.Call(Dx, -layer);
        pc.Call(Dx,  layer);
      }
      for (int Dy = -layer + 1; Dy <= layer - 1; Dy++) {
        pc.Call(-layer, Dy);
        pc.Call( layer, Dy);
      }
    } else if (pE) {
      pc.Call(-layer, 0);
      pc.Call( layer, 0);
    } else {
      pc.Call(0, -layer);
      pc.Call(0,  layer);
    }
  }
  return pc.B();
}

int BrickMeshBemBuilder::SetOmpNthreads (size_t n) {
  omp_set_num_threads(n);
  n = omp_get_max_threads();
  assert (n >= 1);
  if (n > _w.size()) _w.resize(n);
  return n;
}

double BrickMeshBemBuilder::Call (RectId rcv, Dir::Enum dir) {
  assert (Dir::IsCardinal(dir));
  if (_ma.GetBoundaries().GetBC(Dir::ToCardinal(dir))
      != Boundaries::bc_velocity) {
    // There is no point in calling this method for anything other than
    // bc_velocity, but there's no harm.
    return 0;
  }
  assert (omp_get_thread_num() < (int) _w.size());
  return Call_const(rcv, _nr + (int) dir, 0);
}

void BrickMeshBemBuilder::StartBlock () {
  int tid = omp_get_thread_num();
  _w[tid].memo.Clear();
}

void BrickMeshBemBuilder::FinishBlock () {
  int tid = omp_get_thread_num();
//#pragma omp critical (nrgfs_update)
  _ctr.nrgfs += _w[tid].memo.GetNrgfs();
}

inline bool BrickMeshBemBuilder::RcvSrcKey::
operator< (const RcvSrcKey& k) const {
  if (rx < k.rx) return true;
  if (rx > k.rx) return false;
  if (ry < k.ry) return true;
  if (ry > k.ry) return false;
  if (sx < k.sx) return true;
  if (sx > k.sx) return false;
  if (sy < k.sy) return true;
  if (sy > k.sy) return false;
  return false;
}

inline bool BrickMeshBemBuilder::RcvSrcKey::
operator== (const RcvSrcKey& k) const {
  return rx == k.rx && ry == k.ry && sx == k.sx && sy == k.sy;
}

#ifdef USE_HASH_MAP
inline size_t BrickMeshBemBuilder::Memoizer::RcvSrcHasher::
operator() (const RcvSrcKey& k) const {
  // This gem comes from the boost::hash_combine doc. Thanks for the lift,
  // boost! According to some web commentors, STL made a major mistake by not
  // including hash_combine.
  size_t seed;
  const double* vs = &k.rx;
  for (size_t i = 0; i < 4; i++)
    seed ^= hash<double>()(vs[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}
#endif

inline void BrickMeshBemBuilder::Memoizer::Clear () {
  _nrgfs = 0;
  _map.clear();
  // I profiled doing _map.rehash(ceil(_reserve/_map.max_load_factor())), where
  // _rserve is the max observed size up to now, but that's worse.
}

inline pair<BrickMeshBemBuilder::Memoizer::Map::iterator, bool>
BrickMeshBemBuilder::Memoizer::
Lookup (const Rect& sr, double rx, double ry, double* gf) {
  RcvSrcKey key(sr.x + 0.5*sr.dx, sr.y + 0.5*sr.dy, rx, ry);
  pair<Map::iterator, bool> ret = _map.insert(pair<RcvSrcKey, double>(key, 0));
  if (!ret.second) {
    _nrgfs++;
    *gf = ret.first->second;
  }
  return ret;
}

inline void BrickMeshBemBuilder::Memoizer::
Insert (const pair<Map::iterator, bool>& p, double gf) { p.first->second = gf; }

BrickMeshBemBuilder*
NewBrickMeshBemBuilder (const GreensFn* gf, const MeshAnalyzer* ma,
                        const InterpolatorMatrix* im, double dist_fac,
                        size_t n_per_layers)
{ return new BrickMeshBemBuilder(gf, ma, im, dist_fac, n_per_layers); }

BrickMeshBemBuilder*
NewBrickMeshBemBuilder (const GreensFn* gf, const MeshAnalyzer* ma,
                        const InterpolatorMatrix* im, const string& filename)
  throw (FileException)
{ return new BrickMeshBemBuilder(gf, ma, im, filename); }

void DeleteBrickMeshBemBuilder (BrickMeshBemBuilder* bmb) { delete bmb; }

#undef __incr__
#undef incr_ctr
}}
