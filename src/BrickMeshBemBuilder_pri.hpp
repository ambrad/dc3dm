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
 
#ifndef INCLUDE_UTIL_BRICKMESHBEMBUILDER_PRI
#define INCLUDE_UTIL_BRICKMESHBEMBUILDER_PRI

// I find that map is actually a little faster than unordered_map, though I
// thought the opposite would be true.
//#define USE_HASH_MAP
#ifdef USE_HASH_MAP
#include <unordered_map>
#else
#include <map>
#endif
#include "Rect.hpp"
#include "RectMeshUnstruct_pri.hpp"
#include "MeshAnalyzer.hpp"

namespace util {
namespace rmesh {
using namespace std;

// Provide a Green's function for rectangles to points.
class GreensFn {
public:
  virtual ~GreensFn() {}
  virtual double Call(const Rect& src, double rcv_x, double rcv_y) const = 0;
};

// -----------------------------------------------------------------------------
// Builder.
//   [D1] Periodic boundary conditions are handled in the following
// way. Regardless of neighborhood size, interactions between all periodic
// copies of a source except the one closest to the receiver are treated as
// rsr_Normal. This treatment simplifies the code (it lets us keep our upper
// bounds on set sizes that hold for non-periodic-BC problems). Moreover, in
// practice neighborhood sizes shouldn't come close to large enough to make this
// matter.

// Do breadth-first search.
class RmeshBfSearcher {
public:
  RmeshBfSearcher(const MeshAnalyzer* ma);

  // Find (1) RectIds within a unit-cell distance (*not* dist_fac) of root and
  // (2) a supporting set of RectIds beyond these having layer thickness
  // support_layers.
  void FindRectIdsByDist(RectId root, double dist, size_t support_layers);

  // Find RectIds that correspond to rects within nhops of root. Only
  // GetAll() is relevant after this call.
  void FindRectsByHops(RectId root, size_t nhops);

  // The results are ordered
  //     [neighbors supports].
  // GetRectsIds()[GetSupportsStartIdx()] is the first neighbor id.
  size_t GetSupportsStartIdx () const { return _rs.size() - _supports.size(); }
  const vector<RectId>& GetAll () const { return _rs; }
  const vector<RectId>& GetSupport () const { return _supports; }

private:
  const MeshAnalyzer& _ma;
  vector<RectId> _rs, _supports;
  vector<RectId> _discarded; // To clear _in.
  vector<bool> _in; // Mask of RectId in the lists.

  void Clear();
  int Distance(RectId root_id, RectId id) const;
};

// Order in increasing order of strictness.
enum RcvSrcRelation { rsr_Normal = 0, rsr_InterpOnly, rsr_InNbrhd };

struct RcvNbrhd {
  RectId id; // Receiver id.
  size_t nhsz;
  size_t nc; // factor of unit cell (in each dimension)
  // These cells contribute to interpolation but are not subdivided.
  vector<RectId> ins;
  RcvSrcRelation bdy_relations[4];
  // Are periodic sources always rsr_Normal or do they copy the primary source's
  // relation?
  //todo-periodic This is much more efficient than just always using the primary
  // relation, but it isn't optimal.
  bool periodic_normal;
};

struct SrcNbrhd {
  RectId id;
  vector<RectId> ns; // only rects
  vector<RectId> fbps; // bc_free bdy pts
  vector<double> fbps_wts;
};

struct BdyData {
  RectId id;
  // _sns[bdy_pts[i]] is a point on this boundary.
  vector<size_t> bdy_pts;
  Rect r;
};

class BrickMeshBemBuilder {
public:
  BrickMeshBemBuilder(const GreensFn* gf, const MeshAnalyzer* ma,
                      const InterpolatorMatrix* im, double dist_fac,
                      size_t n_per_layers = 0);
  BrickMeshBemBuilder(const GreensFn* gf, const MeshAnalyzer* ma,
                      const InterpolatorMatrix* im, const string& filename)
    throw (FileException);
  ~BrickMeshBemBuilder();

  // Write the neighborhood data. Booleans and OpenMP state are omitted.
  void Serialize(const string& filename) const throw (FileException);
  void Serialize(FILE* fid) const throw (FileException);

  // Call SetOmpNthreads before Call.
  int SetOmpNthreads(size_t n);
  void UsePoorQuality () { _full_quality = false; }
  void UseFullQuality () { _full_quality = true; }
  void SetUseLra (bool use_lra) { _use_lra = use_lra; }
  // tol is optional. If provided, it tells me that I can compute the GF with
  // an absolute error of tol.
  double Call(RectId rcvr, RectId src, double tol = 0);
  double Call(RectId rcvr, Dir::Enum dir);
  void StartBlock();
  void FinishBlock();

  size_t GetBcPeriodicNlayers () const { return _n_per_layers; }
  // Number of calls to Call.
  size_t GetCtrNcalls() const;
  // Number of Green's function evaluations.
  size_t GetCtrNgfs() const;
  size_t GetCtrNrgfs() const;
  double CallGreensFn(const Rect& sr, double rx, double ry,
                      bool use_memo = true) const;

private:
  struct RcvSrcKey {
    double rx, ry, sx, sy;
    RcvSrcKey(double sx, double sy, double rx, double ry)
      : rx(rx), ry(ry), sx(sx), sy(sy) {}
    bool operator<(const RcvSrcKey& k) const;
    bool operator==(const RcvSrcKey& k) const;
  };

  class Memoizer {
  private:
    size_t _nrgfs;

#ifdef USE_HASH_MAP
    struct RcvSrcHasher {
      size_t operator()(const RcvSrcKey& k) const;
    };
    typedef unordered_map<RcvSrcKey, double, RcvSrcHasher> Map;
#else
    typedef map<RcvSrcKey, double> Map;
#endif
    Map _map;

  public:
    typedef Map::iterator iterator;

    Memoizer () : _nrgfs(0) {}
    void Clear();
    size_t GetNrgfs () const { return _nrgfs; }
    pair<Map::iterator, bool> Lookup(const Rect& sr, double rx, double ry,
                                     double* gf);
    void Insert(const pair<Map::iterator, bool>& p, double gf);
  };

public:
  struct Work {
    vector<double> xi, yi, zi;
    vector<Rect> srcs;
    vector<TriTag> in;
    void ClearSrcs();

    Memoizer memo;
  };

private:
  // In this class, RectIds are allocated as
  //     [interior-rects bdy-pts vBC-rects].
  // The third part is new relative to MeshAnalyzer.
  const GreensFn& _gf;
#ifdef SHEAR_SD_FACTOR
  MeshAnalyzer _ma;
#else
  const MeshAnalyzer& _ma;
#endif
  const InterpolatorMatrix& _im;
  bool _any_periodic;
  size_t _nr, _nx;
  vector<RcvNbrhd> _rns;
  vector<SrcNbrhd> _sns;
  BdyData _bdys[4]; // order according to Dir::Enum
  size_t _n_per_layers;
  bool _full_quality, _use_lra;
  // Do a little extra book keeping, even though it's not algorithmically
  // necessary, to make the uniform-mesh (and dist_fac == 0 in general) case a
  // little faster, use less memory, and have a smaller .bmb file.
  bool _dist_fac_is_0;

  mutable vector<Work> _w; // _w[omp_get_thread_num()]
  mutable struct Counter {
    size_t ngfs, ncalls, nrgfs;
    Counter() : ngfs(0), ncalls(0), nrgfs(0) {}
  } _ctr;

private:
  void Deserialize(FILE* fid);
  void BuildBdyData();
  void CalcFbpsWts(const SrcNbrhd& sn);
  void GetWeightedSources(const RcvNbrhd& rn, RectId src) const;
  void GetWeightedSources(const RcvNbrhd& rn, const SrcNbrhd& sn) const;
  bool GetWeightedSourcesSni(const RcvNbrhd& rn, const SrcNbrhd& sn,
                             const SrcNbrhd& sni, Work& w) const;
  RcvSrcRelation GetRelation(const RcvNbrhd& rn, RectId src_id) const;
  Rect GetPrimarySrcRect(RectId rcv, RectId src, Dir::Enum* dir = NULL) const;
  double CallOnPeriodicImages(const Rect& psr, double rx, double ry,
                              bool use_memo = true) const;
  // I want to assert constness even though I have to implement a non-const API.
  double Call_const(RectId rcvr, RectId src, double tol) const;
  double CallHmmvp(RectId rcv, RectId src, double tol) const;
  friend class HmmvpBlockGreensFn;
};

BrickMeshBemBuilder*
NewBrickMeshBemBuilder(const GreensFn* gf, const MeshAnalyzer* ma,
                       const InterpolatorMatrix* im, double dist_fac,
                       size_t n_per_layers = 0);
BrickMeshBemBuilder*
NewBrickMeshBemBuilder(const GreensFn* gf, const MeshAnalyzer* ma,
                       const InterpolatorMatrix* im, const string& filename)
  throw (FileException);
void DeleteBrickMeshBemBuilder(BrickMeshBemBuilder* bmb);
  
}}

#endif
