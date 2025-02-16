#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <climits>

struct IdxDist {
  int idx;
  double dist;
};

// Hashing voxel coordinates
struct VoxKey {
  int x, y, z;
};

struct VoxKeyHash {
  std::size_t operator()(const VoxKey &k) const {
    // Simple hash combination
    // Large primes and a standard hash combine approach
    // Safe for reasonable integer voxel ranges.
    const std::size_t h1 = std::hash<int>()(k.x);
    const std::size_t h2 = std::hash<int>()(k.y);
    const std::size_t h3 = std::hash<int>()(k.z);
    // Combine them
    std::size_t hash_val = h1;
    hash_val = hash_val * 31 + h2;
    hash_val = hash_val * 31 + h3;
    return hash_val;
  }
};

struct VoxKeyEq {
  bool operator()(const VoxKey &a, const VoxKey &b) const {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  }
};

// [[Rcpp::export]]
Rcpp::List radius_search_3d_nonisotropic(
    Rcpp::IntegerMatrix cds_vox,    // Nx3 voxel coords
    Rcpp::NumericMatrix cds_mm,     // Nx3 mm coords
    Rcpp::NumericMatrix queries_mm, // Mx3 query coords in mm
    double radius_mm,
    double sx, double sy, double sz,
    double ox_mm, double oy_mm, double oz_mm
) {
  if (cds_vox.ncol() != 3) Rcpp::stop("cds_vox must be Nx3.");
  if (cds_mm.ncol() != 3) Rcpp::stop("cds_mm must be Nx3.");
  if (queries_mm.ncol() != 3) Rcpp::stop("queries_mm must be Mx3.");

  int n = cds_vox.nrow();
  int qn = queries_mm.nrow();

  if (radius_mm < 0) Rcpp::stop("radius must be non-negative.");

  double r2 = radius_mm * radius_mm;

  // Build hash map from voxel coords -> vector of indices
  std::unordered_map<VoxKey, std::vector<int>, VoxKeyHash, VoxKeyEq> pointMap;
  pointMap.reserve(n*2);
  pointMap.max_load_factor(0.7);

  for (int i = 0; i < n; i++) {
    VoxKey vk {(int)cds_vox(i,0), (int)cds_vox(i,1), (int)cds_vox(i,2)};
    pointMap[vk].push_back(i+1); // store 1-based index
  }

  // Precompute voxel extents for radius search
  // We'll search a voxel cube that definitely covers the mm radius sphere.
  // Radius in voxel units: along each axis separately
  // We expand search radius by converting mm radius to voxel units:
  // For x dimension: #voxels in half-distance = radius_mm / sx, etc.
  Rcpp::List out_indices(qn);
  Rcpp::List out_distances(qn);

  for (int qi = 0; qi < qn; qi++) {
    double qx_mm = queries_mm(qi,0);
    double qy_mm = queries_mm(qi,1);
    double qz_mm = queries_mm(qi,2);

    // Convert query mm to voxel coordinates (as doubles)
    double qx_vox_d = (qx_mm - ox_mm) / sx;
    double qy_vox_d = (qy_mm - oy_mm) / sy;
    double qz_vox_d = (qz_mm - oz_mm) / sz;

    // Determine bounding box in voxel space
    int rx_vox = (int)std::ceil(radius_mm / sx);
    int ry_vox = (int)std::ceil(radius_mm / sy);
    int rz_vox = (int)std::ceil(radius_mm / sz);

    int qx_vox_c = (int)std::floor(qx_vox_d + 0.5);
    int qy_vox_c = (int)std::floor(qy_vox_d + 0.5);
    int qz_vox_c = (int)std::floor(qz_vox_d + 0.5);

    int minX = qx_vox_c - rx_vox;
    int maxX = qx_vox_c + rx_vox;
    int minY = qy_vox_c - ry_vox;
    int maxY = qy_vox_c + ry_vox;
    int minZ = qz_vox_c - rz_vox;
    int maxZ = qz_vox_c + rz_vox;

    std::vector<IdxDist> found_pairs;
    found_pairs.reserve(128);
    std::unordered_set<int> seen_indices;
    seen_indices.reserve(128);

    // Check all voxels in the bounding box
    for (int xx = minX; xx <= maxX; xx++) {
      for (int yy = minY; yy <= maxY; yy++) {
        for (int zz = minZ; zz <= maxZ; zz++) {
          VoxKey vk {xx, yy, zz};
          auto it = pointMap.find(vk);
          if (it != pointMap.end()) {
            // For each point in this voxel, compute mm distance
            const std::vector<int> &indices = it->second;
            for (int idx : indices) {
              if (seen_indices.find(idx) == seen_indices.end()) {
                seen_indices.insert(idx);
                // get mm coords of this point
                double px_mm = cds_mm(idx-1,0);
                double py_mm = cds_mm(idx-1,1);
                double pz_mm = cds_mm(idx-1,2);

                double dx = px_mm - qx_mm;
                double dy = py_mm - qy_mm;
                double dz = pz_mm - qz_mm;
                double dist2 = dx*dx + dy*dy + dz*dz;
                if (dist2 <= r2) {
                  double dist = std::sqrt(dist2);
                  found_pairs.push_back({idx, dist});
                }
              }
            }
          }
        }
      }
    }

    // Sort found points by distance
    std::sort(found_pairs.begin(), found_pairs.end(), [](const IdxDist &a, const IdxDist &b){
      return a.dist < b.dist;
    });

    // Convert to R vectors
    int fcount = (int)found_pairs.size();
    Rcpp::IntegerVector iv(fcount);
    Rcpp::NumericVector dv(fcount);
    for (int i = 0; i < fcount; i++) {
      iv[i] = found_pairs[i].idx;
      dv[i] = found_pairs[i].dist;
    }

    out_indices[qi] = iv;
    out_distances[qi] = dv;
  }

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("indices") = out_indices,
    Rcpp::Named("distances") = out_distances
  );

  return result;
}



// [[Rcpp::export]]
Rcpp::List radius_search_3d_direct(
      Rcpp::IntegerMatrix cds_vox,    // Nx3 voxel coords (ix, iy, iz)
      Rcpp::NumericMatrix cds_mm,     // Nx3 mm coords (x_mm, y_mm, z_mm)
      Rcpp::NumericMatrix queries_mm, // Mx3 queries in mm
      double radius_mm,
      double sx, double sy, double sz,
      double ox_mm, double oy_mm, double oz_mm
  ) {
    if (cds_vox.ncol() != 3) Rcpp::stop("cds_vox must be Nx3.");
    if (cds_mm.ncol() != 3) Rcpp::stop("cds_mm must be Nx3.");
    if (queries_mm.ncol() != 3) Rcpp::stop("queries_mm must be Mx3.");

    int n = cds_vox.nrow();
    int qn = queries_mm.nrow();

    // Compute bounding box
    int xmin = cds_vox(0,0), ymin = cds_vox(0,1), zmin = cds_vox(0,2);
    int xmax = xmin, ymax = ymin, zmax = zmin;
    for (int i = 1; i < n; i++) {
      int x = cds_vox(i,0);
      int y = cds_vox(i,1);
      int z = cds_vox(i,2);
      if (x < xmin) xmin = x;
      if (y < ymin) ymin = y;
      if (z < zmin) zmin = z;
      if (x > xmax) xmax = x;
      if (y > ymax) ymax = y;
      if (z > zmax) zmax = z;
    }

    int dimX = xmax - xmin + 1;
    int dimY = ymax - ymin + 1;
    int dimZ = zmax - zmin + 1;
    size_t voxelCount = (size_t)dimX * (size_t)dimY * (size_t)dimZ;

    // Create a vector of vectors to hold point indices
    std::vector<std::vector<int>> voxelPoints(voxelCount);

    auto voxelIndex = [&](int x, int y, int z) {
      return (size_t)(x - xmin) + ( (size_t)(y - ymin)*dimX ) + ( (size_t)(z - zmin)*dimX*dimY );
    };

    // Fill voxelPoints
    for (int i = 0; i < n; i++) {
      int x = cds_vox(i,0);
      int y = cds_vox(i,1);
      int z = cds_vox(i,2);
      size_t idx = voxelIndex(x,y,z);
      voxelPoints[idx].push_back(i+1); // 1-based index
    }

    double r2 = radius_mm * radius_mm;

    Rcpp::List out_indices(qn);
    Rcpp::List out_distances(qn);

    struct IdxDist {
      int idx;
      double dist;
    };

    for (int qi = 0; qi < qn; qi++) {
      double qx_mm = queries_mm(qi,0);
      double qy_mm = queries_mm(qi,1);
      double qz_mm = queries_mm(qi,2);

      double qx_vox_d = (qx_mm - ox_mm)/sx;
      double qy_vox_d = (qy_mm - oy_mm)/sy;
      double qz_vox_d = (qz_mm - oz_mm)/sz;

      int qx_vox = (int)std::floor(qx_vox_d + 0.5);
      int qy_vox = (int)std::floor(qy_vox_d + 0.5);
      int qz_vox = (int)std::floor(qz_vox_d + 0.5);

      int rx_vox = (int)std::ceil(radius_mm / sx);
      int ry_vox = (int)std::ceil(radius_mm / sy);
      int rz_vox = (int)std::ceil(radius_mm / sz);

      int minX = qx_vox - rx_vox; if (minX < xmin) minX = xmin;
      int maxX = qx_vox + rx_vox; if (maxX > xmax) maxX = xmax;
      int minY = qy_vox - ry_vox; if (minY < ymin) minY = ymin;
      int maxY = qy_vox + ry_vox; if (maxY > ymax) maxY = ymax;
      int minZ = qz_vox - rz_vox; if (minZ < zmin) minZ = zmin;
      int maxZ = qz_vox + rz_vox; if (maxZ > zmax) maxZ = zmax;

      std::vector<IdxDist> found;
      found.reserve(128);

      // Iterate over candidate voxels
      for (int xx = minX; xx <= maxX; xx++) {
        for (int yy = minY; yy <= maxY; yy++) {
          for (int zz = minZ; zz <= maxZ; zz++) {
            size_t vidx = voxelIndex(xx,yy,zz);
            const auto &indices = voxelPoints[vidx];
            for (int idx : indices) {
              double px_mm = cds_mm(idx-1,0);
              double py_mm = cds_mm(idx-1,1);
              double pz_mm = cds_mm(idx-1,2);
              double dx = px_mm - qx_mm;
              double dy = py_mm - qy_mm;
              double dz = pz_mm - qz_mm;
              double dist2 = dx*dx + dy*dy + dz*dz;
              if (dist2 <= r2) {
                found.push_back({idx, std::sqrt(dist2)});
              }
            }
          }
        }
      }

      // Sort by distance
      std::sort(found.begin(), found.end(), [](const IdxDist &a, const IdxDist &b){
        return a.dist < b.dist;
      });

      Rcpp::IntegerVector iv(found.size());
      Rcpp::NumericVector dv(found.size());
      for (size_t i = 0; i < found.size(); i++) {
        iv[(int)i] = found[i].idx;
        dv[(int)i] = found[i].dist;
      }

      out_indices[qi] = iv;
      out_distances[qi] = dv;
    }

    return Rcpp::List::create(
      Rcpp::Named("indices") = out_indices,
      Rcpp::Named("distances") = out_distances
    );
}


  // Precompute all voxel offsets within the radius.
  // Returns a vector of (dx, dy, dz) offsets for voxels within the radius.
  // Checks for radius and spacing > 0 to avoid nonsensical scenarios.
std::vector<std::array<int,3>> precomputeSphereOffsets(double radius_mm, double sx, double sy, double sz) {
    if (radius_mm < 0 || sx <= 0 || sy <= 0 || sz <= 0) {
      Rcpp::stop("Invalid radius or voxel spacing. They must be positive.");
    }

    int rx = (int)std::ceil(radius_mm / sx);
    int ry = (int)std::ceil(radius_mm / sy);
    int rz = (int)std::ceil(radius_mm / sz);

    double r2 = radius_mm * radius_mm;

    std::vector<std::array<int,3>> offsets;
    offsets.reserve((2*rx+1)*(2*ry+1)*(2*rz+1));
    for (int dx = -rx; dx <= rx; dx++) {
      double dx_mm2 = (dx * sx)*(dx * sx);
      for (int dy = -ry; dy <= ry; dy++) {
        double dy_mm2 = (dy * sy)*(dy * sy);
        double dxy = dx_mm2 + dy_mm2;
        if (dxy > r2) continue;
        for (int dz = -rz; dz <= rz; dz++) {
          double dz_mm2 = (dz * sz)*(dz * sz);
          double dist2 = dxy + dz_mm2;
          if (dist2 <= r2) {
            offsets.push_back({dx, dy, dz});
          }
        }
      }
    }
    return offsets;
  }

// [[Rcpp::export]]
Rcpp::List radius_search_3d_precomputed(
    Rcpp::IntegerMatrix cds_vox,    // Nx3 voxel coords
    Rcpp::NumericMatrix cds_mm,     // Nx3 mm coords
    Rcpp::NumericMatrix queries_mm, // Mx3 query coords in mm
    double radius_mm,
    double sx, double sy, double sz,
    double ox_mm, double oy_mm, double oz_mm
) {
  // Basic input checks
  if (cds_vox.ncol() != 3) {
    Rcpp::stop("cds_vox must be Nx3.");
  }
  if (cds_mm.ncol() != 3) {
    Rcpp::stop("cds_mm must be Nx3.");
  }
  if (queries_mm.ncol() != 3) {
    Rcpp::stop("queries_mm must be Mx3.");
  }
  if (cds_vox.nrow() != cds_mm.nrow()) {
    Rcpp::stop("cds_vox and cds_mm must have the same number of rows.");
  }
  if (radius_mm < 0 || sx <= 0 || sy <= 0 || sz <= 0) {
    Rcpp::stop("Invalid radius or voxel spacing. They must be positive.");
  }

  int n = cds_vox.nrow();
  int qn = queries_mm.nrow();

  // Compute bounding box of voxel coordinates
  int xmin = cds_vox(0,0), ymin = cds_vox(0,1), zmin = cds_vox(0,2);
  int xmax = xmin, ymax = ymin, zmax = zmin;
  for (int i = 1; i < n; i++) {
    int x = cds_vox(i,0);
    int y = cds_vox(i,1);
    int z = cds_vox(i,2);
    if (x < xmin) xmin = x;
    if (y < ymin) ymin = y;
    if (z < zmin) zmin = z;
    if (x > xmax) xmax = x;
    if (y > ymax) ymax = y;
    if (z > zmax) zmax = z;
  }

  // Compute dimensions of the bounding box
  // Check for integer overflow - though unlikely for typical volumes
  int64_t dimX = (int64_t)(xmax - xmin + 1);
  int64_t dimY = (int64_t)(ymax - ymin + 1);
  int64_t dimZ = (int64_t)(zmax - zmin + 1);

  if (dimX <= 0 || dimY <= 0 || dimZ <= 0) {
    Rcpp::stop("Invalid bounding box dimensions.");
  }

  // Check for possible overflow in voxelCount
  int64_t voxelCount = dimX * dimY * dimZ;
  if (voxelCount > INT_MAX) {
    Rcpp::stop("Volume too large, voxelCount exceeds INT_MAX.");
  }

  std::vector<std::vector<int>> voxelPoints((size_t)voxelCount);

  auto voxelIndex = [&](int x, int y, int z) {
    int64_t ix = x - xmin;
    int64_t iy = y - ymin;
    int64_t iz = z - zmin;
    int64_t idx = ix + iy*dimX + iz*dimX*dimY;
    return (size_t)idx;
  };

  // Populate voxelPoints
  for (int i = 0; i < n; i++) {
    int x = cds_vox(i,0);
    int y = cds_vox(i,1);
    int z = cds_vox(i,2);
    // Since x,y,z are within [xmin..xmax], [ymin..ymax], [zmin..zmax],
    // voxelIndex should not go out of range.
    size_t idx = voxelIndex(x,y,z);
    voxelPoints[idx].push_back(i+1); // store 1-based index
  }

  // Precompute sphere offsets
  std::vector<std::array<int,3>> offsets = precomputeSphereOffsets(radius_mm, sx, sy, sz);
  double r2 = radius_mm * radius_mm;

  Rcpp::List out_indices(qn);
  Rcpp::List out_distances(qn);

  struct IdxDist {
    int idx;
    double dist2;
  };

  for (int qi = 0; qi < qn; qi++) {
    double qx_mm = queries_mm(qi,0);
    double qy_mm = queries_mm(qi,1);
    double qz_mm = queries_mm(qi,2);

    // Convert query mm to voxel coordinates
    double qx_vox_d = (qx_mm - ox_mm)/sx;
    double qy_vox_d = (qy_mm - oy_mm)/sy;
    double qz_vox_d = (qz_mm - oz_mm)/sz;

    int qx_vox = (int)std::floor(qx_vox_d + 0.5);
    int qy_vox = (int)std::floor(qy_vox_d + 0.5);
    int qz_vox = (int)std::floor(qz_vox_d + 0.5);

    std::vector<IdxDist> found;
    found.reserve(128);

    // Check offsets
    for (auto &off : offsets) {
      int xx = qx_vox + off[0];
      int yy = qy_vox + off[1];
      int zz = qz_vox + off[2];
      // Boundary check before indexing voxelPoints
      if (xx < xmin || xx > xmax || yy < ymin || yy > ymax || zz < zmin || zz > zmax) {
        continue;
      }

      size_t vidx = voxelIndex(xx, yy, zz);
      const auto &indices = voxelPoints[vidx];
      for (int idx : indices) {
        // idx is 1-based, must ensure idx-1 < cds_mm.nrow().
        if (idx < 1 || idx > n) {
          Rcpp::stop("Index out of range: idx in voxelPoints is invalid.");
        }
        double px_mm = cds_mm(idx-1,0);
        double py_mm = cds_mm(idx-1,1);
        double pz_mm = cds_mm(idx-1,2);

        double dx = px_mm - qx_mm;
        double dy = py_mm - qy_mm;
        double dz = pz_mm - qz_mm;
        double dist2 = dx*dx + dy*dy + dz*dz;
        if (dist2 <= r2) {
          found.push_back({idx, dist2});
        }
      }
    }

    std::sort(found.begin(), found.end(), [](const IdxDist &a, const IdxDist &b){
      return a.dist2 < b.dist2;
    });

    Rcpp::IntegerVector iv((int)found.size());
    Rcpp::NumericVector dv((int)found.size());
    for (size_t i = 0; i < found.size(); i++) {
      iv[(int)i] = found[i].idx;
      dv[(int)i] = found[i].dist2;
    }

    out_indices[qi] = iv;
    out_distances[qi] = dv;
  }

  return Rcpp::List::create(
    Rcpp::Named("indices") = out_indices,
    Rcpp::Named("distances") = out_distances
  );
}
