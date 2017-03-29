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
 
#ifndef INCLUDE_UTIL_RECTMESHUNSTRUCT_PRI
#define INCLUDE_UTIL_RECTMESHUNSTRUCT_PRI

#include <vector>
#include "util/include/Matrix.hpp"
#include "util/include/Exception.hpp"
#include "Rect.hpp"
#include "Tri2.hpp"

namespace util {
namespace rmesh {
using namespace std;

// ---------------------------------------------------------------------------
// Mesh the rectangular domain based on the resolution function f(x,y). We
// call the result an f-satisfying mesh.

struct RectOpts {
  // Min and max side lengths.
  double min_len, max_len;
  // Make a uniform mesh with rectangles sized < max_len and respecting
  // max_aspect_ratio.
  bool uniform;
  // Max aspect ratio of a cell.
  static const double max_aspect_ratio;

  RectOpts() : min_len(0.0), max_len(0.0), uniform(false) {}
  RectOpts(FILE* fid);
  RectOpts(double imin_len, double imax_len, bool iuniform)
    : min_len(imin_len), max_len(imax_len), uniform(iuniform) {}
  void Serialize(FILE* fid) const;
  void Deserialize(FILE* fid);
};

class ResolutionFn {
public:
  virtual ~ResolutionFn() {}
  // Get min_{x,y} f for (x,y) \in domain.
  virtual double GetMin(const Rect& domain) const = 0;
};

class TensorMeshLinInterpRF : public ResolutionFn {
public:
  TensorMeshLinInterpRF(const Matrix<double>& x, const Matrix<double>& y,
                        const Matrix<double>& f);

  virtual double GetMin(const Rect& domain) const;

  void GetMin(double* x, double* y, double* f) const;
  void Call(const vector<double>& x, const vector<double>& y,
            vector<double>& f);
  const Matrix<double>& GetX() const;
  const Matrix<double>& GetY() const;
  const Matrix<double>& GetF() const;

private:
  Matrix<double> _x, _y, _f;
};

// Mesh a rectangular domain based on a function f(x, y).
class RectMeshUnstruct {
  friend class RmuAnalyzer;
  friend class RmuSmoother;

public:
  RectMeshUnstruct(const Rect& domain, const RectOpts& ro,
                   ResolutionFn* rf);
  RectMeshUnstruct(const string& filename) throw (FileException);
  ~RectMeshUnstruct();

  void Serialize(const string& filename) const throw (FileException);
  void Serialize(FILE* fid) const throw (FileException);

  const Rect& GetDomain() const;
  const RectOpts& GetRectOpts() const { return _ro; }
  const vector<Rect>& GetRects() const;
  // Return base-0 index i such that GetRects()[i] is the rectangle that
  // contains the point (x, y). Result is -1 if (x, y) is outside of the
  // domain.
  int GetRectId(double x, double y) const;

  // For convergence testing, refine each element by one level.
  void Refine();

private:
  class QuadTree;
      
  Rect _domain;
  RectOpts _ro;
  vector<QuadTree*> _qts;
  size_t _nx, _ny;
  mutable vector<Rect> _rs;

private:
  void Deserialize(FILE* fid);
  RectMeshUnstruct(const RectMeshUnstruct& rmu);
};

class RectMeshUnstruct::QuadTree {
  friend class RmuAnalyzer;
  friend class RmuSmoother;

private:
  Rect _r;
  union {
    struct {
      QuadTree* kids[4]; // Ordered by quadrant (NE, NW, SW, SE)
    } _n; // Internal node
    struct {
      RectId id;
    } _l; // Leaf
  };
  bool _is_leaf;

public:
  QuadTree(const Rect& r);
  QuadTree(FILE* fid);
  QuadTree(const QuadTree& qt);
  ~QuadTree();

  const Rect& GetRect() const { return _r; }
  void Split();
  void Split(const RectOpts& ro, ResolutionFn* rf, RectId* starting_id,
             bool first = false);
  void PushBackRects(vector<Rect>* rs);
  int GetRectId(double x, double y) const;
  void Refine(RectId* id);
  void Serialize(FILE* fid);
  void Deserialize(FILE* fid);
};

RectMeshUnstruct*
NewRectMeshUnstruct(const Rect& domain, const RectOpts& ro, ResolutionFn* rf);
RectMeshUnstruct* NewRectMeshUnstruct(const string& filename)
  throw (FileException);
RectMeshUnstruct*
SmoothRectMeshUnstruct(const RectMeshUnstruct* rmu, const Boundaries& b);
void DeleteRectMeshUnstruct(RectMeshUnstruct* rmu);

// -----------------------------------------------------------------------------
// Take an f-satisfying mesh and (1) refine it based on boundaries and (2)
// obtain more complex data structures.

struct TriTag {
  size_t id; // Index into tri.
  // Translate at the domain level.
  //todo-periodic-complete dir -> (Dx, Dy)
  Dir::Enum dir;
  // The vertex we want anchored at a rect center. This is used only if there is
  // a periodic BC.
  char anchor;
};

//todo-periodic Unused.
struct RectTag {
  RectId id;
  // Domain shift. For example, if this periodic rect is two domains N and one
  // W of the domain, then Dx = -1 and Dy = 2. Only used if there is a
  // periodic BC.
  char Dx, Dy;
};

struct Edge {
  RectId id;     // Connect to this node.
  Dir::Enum bdy; // Possibly go through a periodic boundary to do so.
  bool use_in_GetEdge; // Handle an issue with v-BC - periodic-BC corners.

  Edge() : id(0), bdy(Dir::inside), use_in_GetEdge(true) {}
  Edge(RectId iid) : id(iid), bdy(Dir::inside), use_in_GetEdge(true) {}
  Edge(RectId iid, Dir::Enum ibdy)
    : id(iid), bdy(ibdy), use_in_GetEdge(true) {}

  // NB: Relations do not use bdy.
  bool operator< ( const Edge& e2) const { return this->id <  e2.id; }
  bool operator== (const Edge& e2) const { return this->id == e2.id; }
};

typedef vector< vector<Edge> > Vve;

class RmuAnalyzer {
public:
  RmuAnalyzer(const RectMeshUnstruct* rmu, const Boundaries& b);
  RmuAnalyzer(const RectMeshUnstruct* rmu, const string& filename)
    throw (FileException);
  ~RmuAnalyzer();

  void Serialize(const string& filename) const throw (FileException);
  void Serialize(FILE* fid) const throw (FileException);

  const Boundaries& GetBoundaries() const;
  // Get the neighbors of rectangle id. A neighbor is either (1) a rectangle,
  // either directly adjacent or periodically adjacent, or (2) a point on a
  // boundary corresponding to a velocity boundary condition.
  const vector<RectId>& GetNbrs(RectId id) const;
  const vector<Dir::Enum>& GetNbrsBdys(RectId id) const;
  // Returns Dir::inside if not on boundary. Only boundary points are considered
  // to be on boundaries in this method; rectangles are not.
  Dir::Enum GetBoundaryDir(RectId id) const;
  // Accounts for periodicity.
  //todo-periodic Not sure this is needed.
  void Distance(RectId id1, RectId id2, double* dx, double* dy,
                // Optionally ask for the boundary one goes through to get
                // from id1 to id2.
                Dir::Enum* bdy = NULL) const;
  void Distance(const RectTag& rt1, const RectTag& rt2, double* dx,
                double* dy) const;
  const Matrix<double>& GetX() const;
  const Matrix<RectId>& GetTri() const;
  // For base-1 index tag->id, tri = GetTri(), tri(:, tag->id) is the triangle
  // containing (x, y). Returns false if (x, y) is not in the domain.
  bool GetTriTag(const double x, const double y, TriTag* tag) const;

  // Get zi[i](xi[i], yi[i]). zi[i] is set to 0 if (xi[i], yi[i]) is outside
  // of the domain of interpolatable points.
  //   First call this method. On entry, z_int[i] is z(center of GetX()[i])
  // for i = 0 to GetRects().size() - 1, ie, the values at the interior
  // points. On exit, it contains boundary information and is suitable for
  // input to LinearInterp.
  //   Provide boundary values for each side; only those on v-BC sides are
  // actually used. Order: E N W S.
  void PrepInterp(vector<double>& z_int, const double bdy_vals[4])
    const;
  //   Alternatively, call this routine to use extrapolation on all but
  // (obviously) periodic boundaries.
  void PrepExtrap(vector<double>& z_int) const;
  void LinearInterp(
    // z[i] is z(center of GetX()[i]) for i = 0 to GetX().Size(2) - 1,
    // probably from PrepLinearInterp.
    const vector<double>& z,
    // Interpolate to these points.
    const vector<double>& xi, const vector<double>& yi,
    // Interpolated values. On output, zi.size() == xi.size().
    vector<double>& zi) const;
  void GetTrisSharingRect(RectId id, vector<size_t>& tis) const;

  // Public in the _pri interface.

  //   [c1] A corner point at the corner of two v-BCs is assigned to one
  // boundary or the other. If the v-BC slabs are slipping the same, this
  // doesn't matter; in my view, it's not correct for them not to be.
  //   A corner point at the corner of a (0-)v-BC and a free BC is assigned to
  // the (0-v)-BC side.
  //   A corner point at the corner of a 0-v-BC and a v-BC is assigned to the
  // v-BC side.
  void GetCornerBC(Dir::Enum corner, Dir::Enum& bdy, Boundaries::BC& bc) const;
  Dir::Enum OnWhichBoundary(const Rect& r) const;
  Dir::Enum OnWhichBoundary(const double x, const double y) const;
  // Periodic mapping to implement boundary conditions. "Periodic" is used to
  // mean the mapping of one point to another based on periodic boundary
  // conditions. "Domain" is the rectangle GetDomain(). int arguments indicate
  // calculations are in unit coordinates rather than domain
  // coordinates. Eventually, I want to do most calcs in unit coords.
  bool MapToDomain(double* x, double* y, Dir::Enum* bdy) const;
  bool MapToCover(double x, double y, Rect* r) const;
  void MapToPeriodicWithTriAnchoredAt(
    size_t anchor_id, size_t v_id, const RectId tri[3],
    const Dir::Enum tri_bdy[3], double* x, double* y) const;
  void MapToPeriodic(Dir::Enum bdy, double* x, double* y) const;
  void MapToPeriodic(Dir::Enum bdy, int* x, int* y) const;
  // Get triangle vertices. The purpose of this routine is to map vertices as
  // necessary for periodic boundary conditions.
  void GetPeriodicTri2(const TriTag& tag, double* v1, double* v2, double* v3)
    const;
  int UnitDx(double dx) const;
  int UnitDy(double dy) const;
  int UnitX(double x) const;
  int UnitY(double y) const;
  double DomainX(int x) const;
  double DomainY(int y) const;

  //testing
  const Matrix<Dir::Enum>& GetTriBdy() const;
  void InterpWExtrap(const vector<double>& z_int, const vector<double>& xi,
                     const vector<double>& yi, vector<double>& zi) const;

private:
  const RectMeshUnstruct& _rmu;
  const Rect& _r;  // The domain.
  const vector<Rect>& _rs;

  Boundaries _b;   // Boundary conditions.
  double _ux, _uy; // Unit x and y.
  size_t _nr;
  bool _any_periodic;

  // Adjacency graph: _ag[id] is a sorted (by increasing id) list of nbrs. ids
  // are divided into two sets: rectangles and boundary points. ids are
  // allocated first to rectangles. RectIds start at 0.
  vector< vector<RectId> > _ag;
  // _ag_bdy[i][j] is the boundary(s) to go through to get from node i, which is
  // assumed to be in the domain, to j, which may be in a periodic rectangle.
  // This relationship is just a convention; the reverse relationship equally
  // holds. If no boundary is traversed, it's set to Dir::inside. If two
  // periodic boundaries are traversed, the direction is noncardinal.
  vector< vector<Dir::Enum> > _ag_bdy;
  // Rect centers + boundary rect corners (for v-BC).
  Matrix<double> _x;
  // _bdypt_side[i - _nr] is the side or corner the edge point is on.
  vector<Dir::Enum> _bdypt_side;
  // Triangulation.
  Matrix<RectId> _tri; // _tri(:,i) is the i'th triangle's RectIds.
  // _tri_bdy(i,j) for vertex i of triangle j says whether the vertex is reached
  // through a boundary. A triangle that crosses a periodic boundary is one
  // particular instance of up to three such triangles. TriTag is used to
  // distinguish among the periodic images.
  Matrix<Dir::Enum> _tri_bdy;
  // For fast triangle finding.
  struct RtoT {
    size_t ti;
    double theta;
    bool operator<(const RtoT& r1) const;
  };
  vector< vector<RtoT> > _r2t;

private:
  void Deserialize(FILE* fid);
  void OrganizeInSpace(Vve& age);
  void OrganizeInSpaceHandleBcCornerVP(
    RectId rect_id, RectId node_id, Dir::Enum bdy, vector<bool>& handled_corner,
    Vve& age);
  void GetNbrs(RectId id, vector<Edge>& nbrs) const;
  void RefineAtBoundaries(double dist_fac);
  void Triangulate(const Vve& age, Matrix<Edge>& trie);
  void MakeRtoT();
  double Angle(RectId ir, const RectId tri[3], const Dir::Enum tri_bdy[3])
    const;
  void GetTriTags(RectId rid, double x, double y, TriTag tag[4]) const;
  bool IsPointInTri(const TriTag& tag, double x, double y) const;
  int IsPointInTris(const TriTag tags[3], double ox, double oy) const;
  double Extrap(const vector<double>& z_int, RectId id) const;

  // For internal use by RmuSmoother, this constructs just _ag, _ag_bdy,
  // _bdypt_side.
  RmuAnalyzer(const RectMeshUnstruct* rmu, const Boundaries& b,
              bool want_ag_only);
  friend RectMeshUnstruct* SmoothRectMeshUnstruct(const RectMeshUnstruct* rmu,
                                                  const Boundaries& b);
};

RmuAnalyzer* NewRmuAnalyzer(const RectMeshUnstruct* rmu, const Boundaries& b);
RmuAnalyzer* NewRmuAnalyzer(const RectMeshUnstruct* rmu,
                            const string& filename)
  throw (FileException);
void DeleteRmuAnalyzer(RmuAnalyzer* ra);

}}

#endif
