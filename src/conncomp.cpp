// Connected-component labelling for 3-D masks.
//
// This replaces a two-pass union-find written in interpreted R, which built a
// 26x3 neighbour matrix per voxel and ran find() as an R closure. Measured
// against scipy.ndimage.label() the R version was 690x slower on a 40x48x38
// volume and 4,450x slower on 182x218x182 -- fifteen minutes to label a 1 mm
// whole brain -- and the ratio grew with size. automask() spent 92% of its
// runtime here, so this is its cost too.
//
// The algorithm and, importantly, the *output numbering* are the same as the R
// version's, so the change is invisible:
//
//   * Voxels are visited in column-major order. Each gets a provisional label
//     equal to its ordinal in that order, whether or not it ends up using it --
//     the R loop incremented `nextlabel` on every voxel.
//   * A voxel with no already-labelled neighbour starts a new set; otherwise it
//     takes the minimum neighbour label ML, and every distinct neighbour label
//     is unioned into find(ML).
//   * Components are then numbered by size descending, ties broken by the
//     smaller provisional root label -- which is what
//     `sort(table(labs), decreasing = TRUE)` produced, since table() orders its
//     names numerically and R's sort is stable.
//
// tests/testthat/test-conncomp-equivalence.R pins the result against the
// reference R implementation kept there, over randomised masks and all three
// connectivities.

#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rcpp.h>
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif

#include <algorithm>
#include <numeric>
#include <vector>

using namespace Rcpp;

namespace {

// Neighbour offsets, in the row order the R implementation used. The order is
// not load-bearing -- the label a voxel takes is the minimum over neighbours,
// and the root a union lands on is find(ML) either way -- but keeping it makes
// the two implementations comparable step by step.
const int OFF6[6][3] = {
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};

const int OFF18[18][3] = {
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1},
    {-1, -1, 0}, {-1, 1, 0}, {1, -1, 0}, {1, 1, 0},
    {-1, 0, -1}, {-1, 0, 1}, {1, 0, -1}, {1, 0, 1},
    {0, -1, -1}, {0, -1, 1}, {0, 1, -1}, {0, 1, 1}};

const int OFF26[26][3] = {
    {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1},
    {-1, 1, -1}, {-1, 1, 0}, {-1, 1, 1},
    {0, -1, -1}, {0, -1, 0}, {0, -1, 1}, {0, 0, -1}, {0, 0, 1},
    {0, 1, -1}, {0, 1, 0}, {0, 1, 1},
    {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, {1, 0, -1}, {1, 0, 0}, {1, 0, 1},
    {1, 1, -1}, {1, 1, 0}, {1, 1, 1}};

struct DisjointSet {
  std::vector<int> parent;

  explicit DisjointSet(size_t n) : parent(n) {
    std::iota(parent.begin(), parent.end(), 0);
  }
  int find(int i) {
    int root = i;
    while (parent[root] != root) root = parent[root];
    while (parent[i] != root) {          // path compression
      int next = parent[i];
      parent[i] = root;
      i = next;
    }
    return root;
  }
  void unite(int a, int b) {
    int ra = find(a), rb = find(b);
    if (ra != rb) parent[ra] = rb;
  }
};

}  // namespace

//' Label the connected components of a 3-D logical mask
//'
//' @param mask logical vector holding the mask in column-major order
//' @param dims integer vector of length 3
//' @param connectivity 6, 18 or 26
//' @return a list with \code{index} (components numbered by decreasing size)
//'   and \code{size} (the size of the component each voxel belongs to), both
//'   integer vectors of \code{prod(dims)}, and \code{n} the component count
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List conn_comp_labels_cpp(LogicalVector mask, IntegerVector dims, int connectivity) {
  if (dims.size() != 3) stop("conn_comp_labels_cpp: 'dims' must have length 3");
  const int d0 = dims[0], d1 = dims[1], d2 = dims[2];
  if (d0 < 0 || d1 < 0 || d2 < 0) stop("conn_comp_labels_cpp: 'dims' must be non-negative");

  const double dnvox = (double) d0 * (double) d1 * (double) d2;
  if (dnvox > (double) R_XLEN_T_MAX) stop("conn_comp_labels_cpp: volume too large");
  const R_xlen_t nvox = (R_xlen_t) dnvox;
  if (mask.size() != nvox)
    stop("conn_comp_labels_cpp: 'mask' has %lld elements but 'dims' implies %lld",
         (long long) mask.size(), (long long) nvox);

  const int (*off)[3];
  int n_off;
  switch (connectivity) {
    case 6:  off = OFF6;  n_off = 6;  break;
    case 18: off = OFF18; n_off = 18; break;
    case 26: off = OFF26; n_off = 26; break;
    default: stop("conn_comp_labels_cpp: 'connectivity' must be 6, 18 or 26");
  }

  IntegerVector index(nvox), size_out(nvox);
  if (nvox == 0) return List::create(_["index"] = index, _["size"] = size_out, _["n"] = 0);

  const int* mp = LOGICAL(mask);
  // Provisional labels are 1-based ordinals over the in-mask voxels, so the
  // set can be no larger than the number of them.
  R_xlen_t n_in = 0;
  for (R_xlen_t i = 0; i < nvox; i++) if (mp[i] == TRUE) n_in++;
  if (n_in == 0) return List::create(_["index"] = index, _["size"] = size_out, _["n"] = 0);
  if (n_in > (R_xlen_t) INT_MAX)
    stop("conn_comp_labels_cpp: more than 2^31 mask voxels is not supported");

  std::vector<int> labels((size_t) nvox, 0);
  DisjointSet ds((size_t) n_in + 1);

  const R_xlen_t slice = (R_xlen_t) d0 * d1;
  std::vector<int> nabe;                       // neighbour labels, in offset order
  nabe.reserve((size_t) n_off);
  int next_label = 0;

  for (int z = 0; z < d2; z++) {
    for (int y = 0; y < d1; y++) {
      for (int x = 0; x < d0; x++) {
        const R_xlen_t v = (R_xlen_t) x + (R_xlen_t) y * d0 + (R_xlen_t) z * slice;
        if (mp[v] != TRUE) continue;
        ++next_label;                          // one per in-mask voxel, as in R

        nabe.clear();
        for (int k = 0; k < n_off; k++) {
          const int nx = x + off[k][0], ny = y + off[k][1], nz = z + off[k][2];
          if (nx < 0 || nx >= d0 || ny < 0 || ny >= d1 || nz < 0 || nz >= d2) continue;
          const R_xlen_t u = (R_xlen_t) nx + (R_xlen_t) ny * d0 + (R_xlen_t) nz * slice;
          const int lab = labels[(size_t) u];
          if (lab != 0) nabe.push_back(lab);
        }

        if (nabe.empty()) {
          labels[(size_t) v] = next_label;     // starts its own set
        } else {
          const int ml = *std::min_element(nabe.begin(), nabe.end());
          labels[(size_t) v] = ml;
          ds.parent[(size_t) next_label] = ml; // the unused ordinal follows ML
          for (size_t j = 0; j < nabe.size(); j++) ds.unite(nabe[j], ml);
        }
      }
    }
  }

  // Resolve every voxel to its root, and count the roots.
  std::vector<int> root_size((size_t) n_in + 1, 0);
  for (R_xlen_t v = 0; v < nvox; v++) {
    if (labels[(size_t) v] == 0) continue;
    const int r = ds.find(labels[(size_t) v]);
    labels[(size_t) v] = r;
    root_size[(size_t) r]++;
  }

  // Number the components: size descending, then root label ascending. That is
  // what sort(table(labs), decreasing = TRUE) produced in R.
  std::vector<int> roots;
  for (int r = 1; r <= (int) n_in; r++) if (root_size[(size_t) r] > 0) roots.push_back(r);
  std::stable_sort(roots.begin(), roots.end(),
                   [&](int a, int b) { return root_size[(size_t) a] > root_size[(size_t) b]; });

  std::vector<int> rank((size_t) n_in + 1, 0);
  for (size_t i = 0; i < roots.size(); i++) rank[(size_t) roots[i]] = (int) i + 1;

  int* ip = INTEGER(index);
  int* sp = INTEGER(size_out);
  for (R_xlen_t v = 0; v < nvox; v++) {
    const int r = labels[(size_t) v];
    if (r == 0) { ip[v] = 0; sp[v] = 0; continue; }
    ip[v] = rank[(size_t) r];
    sp[v] = root_size[(size_t) r];
  }

  return List::create(_["index"] = index, _["size"] = size_out,
                      _["n"] = (int) roots.size());
}

// ---------------------------------------------------------------------------
// Local-maxima pruning
//
// conn_comp(local_maxima = TRUE) thins each component down to peaks that are at
// least `mindist` apart. A point is kept when no other point within mindist has
// a larger value; that is the property `local_maxima_dist` names, and it holds
// after a single pass, because a kept point is defined against the whole input
// rather than against whatever survived alongside it.
//
// What this replaces did something weaker and did it in R: it compared each
// point only against its single nearest neighbour, via a dbscan::kNN() call per
// component per round. Two things were wrong with that. It does not deliver the
// minimum distance it promises -- with mindist = 15, points at 0, 5 and 12 mm
// valued 5, 1 and 9 keep the 5 and the 9, twelve millimetres apart, because the
// 5's *nearest* neighbour is the 1 -- and "the nearest neighbour" is not even
// well defined on a voxel grid, where more than half of all points have several
// neighbours at exactly the same distance and which one came back depended on
// the kd-tree's traversal order.
//
// Cost: points are bucketed into a uniform grid of cells `mindist` wide, so
// anything within mindist of a point lies in one of the 27 cells around it.
// Each point scans those cells and stops at the first larger-valued point close
// enough to knock it out, which is immediate for the many points that are not
// peaks; the full scan is paid only by the few that survive.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector prune_local_maxima_cpp(NumericMatrix coords, NumericVector vals,
                                     double mindist) {
  const int n = coords.nrow();
  const int ndim = coords.ncol();
  if (vals.size() != n) stop("prune_local_maxima_cpp: 'vals' must have nrow(coords) entries");
  if (ISNA(mindist)) stop("prune_local_maxima_cpp: 'mindist' must not be NA");

  // Nothing can be strictly closer than a non-positive distance, so every point
  // is a maximum by default.
  if (n <= 1 || !(mindist > 0.0)) {
    IntegerVector all(n);
    for (int i = 0; i < n; i++) all[i] = i + 1;
    return all;
  }

  // Cell coordinates, one per axis, in cells of side `mindist`. Origin at the
  // bounding-box corner so the indices are non-negative.
  std::vector<double> lo((size_t) ndim, R_PosInf);
  for (int c = 0; c < ndim; c++)
    for (int i = 0; i < n; i++)
      if (coords(i, c) < lo[c]) lo[c] = coords(i, c);

  std::vector<int> cell((size_t) n * ndim);
  std::vector<int> ncell((size_t) ndim, 1);
  for (int c = 0; c < ndim; c++) {
    double hi = 0.0;
    for (int i = 0; i < n; i++) {
      const double t = (coords(i, c) - lo[c]) / mindist;
      if (!R_finite(t)) stop("prune_local_maxima_cpp: 'coords' must be finite");
      cell[(size_t) i * ndim + c] = (int) t;
      if (t > hi) hi = t;
    }
    ncell[c] = (int) hi + 1;
  }

  // A dense cell array would be (span/mindist)^ndim, which a long thin
  // component makes far larger than the point count, so the cells are hashed
  // into n buckets and held as a CSR-style list: `start` indexes `order`.
  const size_t nbuck = (size_t) n;
  std::vector<int> bucket((size_t) n);
  auto hash_cell = [&](const int* cc) {
    unsigned long long h = 1469598103934665603ULL;
    for (int c = 0; c < ndim; c++) {
      h ^= (unsigned long long)(unsigned int) cc[c];
      h *= 1099511628211ULL;
    }
    return (size_t)(h % nbuck);
  };
  std::vector<int> start(nbuck + 1, 0);
  for (int i = 0; i < n; i++) {
    bucket[i] = (int) hash_cell(&cell[(size_t) i * ndim]);
    start[(size_t) bucket[i] + 1]++;
  }
  for (size_t b = 0; b < nbuck; b++) start[b + 1] += start[b];
  std::vector<int> order((size_t) n);
  {
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (int i = 0; i < n; i++) order[(size_t) fill[(size_t) bucket[i]]++] = i;
  }

  // Offsets over the 3^ndim cells around a point.
  std::vector<int> shifts;
  {
    std::vector<int> cur((size_t) ndim, -1);
    for (;;) {
      shifts.insert(shifts.end(), cur.begin(), cur.end());
      int c = ndim - 1;
      while (c >= 0 && cur[(size_t) c] == 1) { cur[(size_t) c] = -1; c--; }
      if (c < 0) break;
      cur[(size_t) c]++;
    }
  }
  const size_t nshift = shifts.size() / (size_t) ndim;

  const double md2 = mindist * mindist;
  std::vector<int> keep;
  std::vector<int> probe((size_t) ndim);

  for (int i = 0; i < n; i++) {
    const double vi = vals[i];
    const int* ci = &cell[(size_t) i * ndim];
    bool beaten = false;

    for (size_t s = 0; s < nshift && !beaten; s++) {
      bool inside = true;
      for (int c = 0; c < ndim; c++) {
        const int p = ci[c] + shifts[s * (size_t) ndim + (size_t) c];
        if (p < 0 || p >= ncell[(size_t) c]) { inside = false; break; }
        probe[(size_t) c] = p;
      }
      if (!inside) continue;

      const size_t b = hash_cell(probe.data());
      for (int k = start[b]; k < start[b + 1] && !beaten; k++) {
        const int j = order[(size_t) k];
        if (j == i || !(vals[j] > vi)) continue;
        // Buckets collide, so confirm the cell before trusting the bucket.
        bool same = true;
        for (int c = 0; c < ndim; c++)
          if (cell[(size_t) j * ndim + c] != probe[(size_t) c]) { same = false; break; }
        if (!same) continue;
        double d2 = 0.0;
        for (int c = 0; c < ndim; c++) {
          const double dd = coords(i, c) - coords(j, c);
          d2 += dd * dd;
        }
        if (d2 < md2) beaten = true;
      }
    }

    if (!beaten) keep.push_back(i + 1);
  }

  return IntegerVector(keep.begin(), keep.end());
}
