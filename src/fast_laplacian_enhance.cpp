#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <RcppArmadillo.h>
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Helper function to mirror index
inline int mirrorIndex(int i, int max_i) {
  if (i < 0) return -i - 1;
  if (i >= max_i) return 2*max_i - i - 1;
  return i;
}

// Get voxel value with mirroring
inline double getVolValue(const NumericVector &vol, int x, int y, int z,
                          int nx, int ny, int nz) {
  x = mirrorIndex(x, nx);
  y = mirrorIndex(y, ny);
  z = mirrorIndex(z, nz);
  return vol[x + nx*(y + ny*z)];
}

// Extract a patch on-demand
static void extract_patch(const NumericVector &vol, int nx, int ny, int nz,
                          int x, int y, int z, int half_p,
                          std::vector<double> &patch_out) {
  int ps = 2*half_p + 1;
  patch_out.resize(ps*ps*ps);
  int idx = 0;
  for (int zz = -half_p; zz <= half_p; zz++) {
    for (int yy = -half_p; yy <= half_p; yy++) {
      for (int xx = -half_p; xx <= half_p; xx++) {
        double val = getVolValue(vol, x+xx, y+yy, z+zz, nx, ny, nz);
        patch_out[idx++] = val;
      }
    }
  }
}

// Build W from K and d (no in-place modification)
static arma::sp_mat buildW(const arma::sp_mat &K, const arma::vec &d_vals, bool use_nf) {
  int n = (int)K.n_rows;
  if (use_nf) {
    double alpha = 1.0/arma::mean(d_vals);

    // alpha*(K - D) + I
    // We'll reconstruct from scratch
    std::vector<arma::uword> Ii;
    std::vector<arma::uword> Jj;
    std::vector<double> Vv;
    Ii.reserve(K.n_nonzero + (size_t)n);
    Jj.reserve(K.n_nonzero + (size_t)n);
    Vv.reserve(K.n_nonzero + (size_t)n);

    // Store entries of alpha*(K-D)
    // K-D means subtract d_vals[r] from diagonal elements.
    // Then add I at the end.
    std::unordered_map<std::size_t,double> diag_map;
    diag_map.reserve((size_t)n);

    for (auto it = K.begin(); it!=K.end(); ++it) {
      arma::uword r = it.row();
      arma::uword c = it.col();
      double val = (*it);
      if (r == c) {
        val -= d_vals[r];
      }
      val *= alpha;
      if (std::fabs(val) > 1e-15) {
        Ii.push_back(r);
        Jj.push_back(c);
        Vv.push_back(val);
        if (r == c) {
          std::size_t key = (std::size_t)r + (std::size_t)n*(std::size_t)c;
          diag_map[key] = (double)Vv.size()-1; // store index
        }
      }
    }

    // Add I
    for (int i=0; i<n; i++) {
      std::size_t key = (std::size_t)i + (std::size_t)n*(std::size_t)i;
      // Check if diagonal was present
      //double added = 1.0;
      // If diagonal element existed in Vv:
      auto f = diag_map.find(key);
      if (f != diag_map.end()) {
        // increment existing diagonal
        arma::uword idx = (arma::uword)f->second;
        Vv[idx] += 1.0;
      } else {
        // add new diagonal
        Ii.push_back((arma::uword)i);
        Jj.push_back((arma::uword)i);
        Vv.push_back(1.0);
      }
    }

    arma::umat locs(2, Ii.size());
    for (size_t i=0; i<Ii.size(); i++){
      locs(0,i) = Ii[i];
      locs(1,i) = Jj[i];
    }
    arma::vec vals(Vv);
    arma::sp_mat W(locs, vals, (arma::uword)n, (arma::uword)n);
    return W;
  } else {
    // W = D^-1 K
    std::vector<arma::uword> Ii;
    std::vector<arma::uword> Jj;
    std::vector<double> Vv;
    Ii.reserve(K.n_nonzero);
    Jj.reserve(K.n_nonzero);
    Vv.reserve(K.n_nonzero);

    for (auto it = K.begin(); it != K.end(); ++it) {
      arma::uword r = it.row();
      arma::uword c = it.col();
      double val = (*it)/(d_vals[r]+1e-15);
      if (std::fabs(val)>1e-15) {
        Ii.push_back(r);
        Jj.push_back(c);
        Vv.push_back(val);
      }
    }
    arma::umat locs(2, Ii.size());
    for (size_t i=0; i<Ii.size(); i++){
      locs(0,i) = Ii[i];
      locs(1,i) = Jj[i];
    }
    arma::vec vals(Vv);
    arma::sp_mat W(locs, vals, (arma::uword)n, (arma::uword)n);
    return W;
  }
}

// Hadamard (element-wise) product for sparse
static arma::sp_mat hadamard(const arma::sp_mat &A, const arma::sp_mat &B) {
  // Create a map for B
  std::unordered_map<std::size_t,double> Bmap;
  Bmap.reserve(B.n_nonzero);
  {
    for (auto it = B.begin(); it != B.end(); ++it) {
      std::size_t key = (std::size_t)it.row() + ((std::size_t)B.n_rows*(std::size_t)it.col());
      Bmap[key] = (*it);
    }
  }

  std::vector<arma::uword> Ii;
  std::vector<arma::uword> Jj;
  std::vector<double> Vv;
  Ii.reserve(std::min(A.n_nonzero, B.n_nonzero));
  Jj.reserve(std::min(A.n_nonzero, B.n_nonzero));
  Vv.reserve(std::min(A.n_nonzero, B.n_nonzero));

  for (auto it = A.begin(); it != A.end(); ++it) {
    arma::uword r = it.row();
    arma::uword c = it.col();
    std::size_t key = (std::size_t)r + ((std::size_t)A.n_rows*(std::size_t)c);
    auto f = Bmap.find(key);
    if (f != Bmap.end()) {
      double val = (*it)*f->second;
      if (std::fabs(val)>1e-15) {
        Ii.push_back(r);
        Jj.push_back(c);
        Vv.push_back(val);
      }
    }
  }

  arma::umat locs(2, Ii.size());
  for (size_t i=0; i<Ii.size(); i++) {
    locs(0,i)=Ii[i]; locs(1,i)=Jj[i];
  }
  arma::vec vals(Vv);
  arma::sp_mat C(locs, vals, A.n_rows, A.n_cols);
  return C;
}

// Sigmoid mappings
inline double detail_map(double x, double a, double width) {
  return 2.0/(1.0+std::exp(-a*x*width))-1.0;
}
inline double base_map(double x, double a, double width) {
  return 1.0/(1.0+std::exp(-a*(x-0.5)*width));
}

// [[Rcpp::export]]
NumericVector fast_multilayer_laplacian_enhancement_masked(
    NumericVector img,
    LogicalVector mask,
    int k=2,
    int patch_size=3,
    int search_radius=2,
    double h=0.7,
    List mapping_params = R_NilValue,
    bool use_normalization_free=true
) {
  if (!img.hasAttribute("dim"))
    stop("img must have dim attribute");
  if (!mask.hasAttribute("dim"))
    stop("mask must have dim attribute");

  IntegerVector dims = img.attr("dim");
  if (dims.size() != 3) stop("img must be 3D");
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];
  int nvox = nx*ny*nz;

  IntegerVector mask_dims = mask.attr("dim");
  if (mask_dims.size()!=3 || mask_dims[0]!=nx || mask_dims[1]!=ny || mask_dims[2]!=nz)
    stop("mask dimensions must match img");

  if (patch_size % 2 == 0)
    stop("patch_size must be odd.");
  if (search_radius <=0)
    stop("search_radius must be positive.");

  double img_min = Rcpp::min(img);
  double img_max = Rcpp::max(img);
  double scale = img_max - img_min + 1e-15;
  NumericVector img_norm(nvox);
  for (int i=0; i<nvox; i++) {
    img_norm[i] = (img[i]-img_min)/scale;
  }

  // Create masked index mapping
  std::vector<int> masked_indices;
  masked_indices.reserve((size_t)nvox);
  for (int i=0; i<nvox; i++) {
    if (mask[i]) {
      masked_indices.push_back(i);
    }
  }
  int M = (int)masked_indices.size();
  if (M==0)
    stop("No voxels inside mask");

  // index_map: maps global idx -> masked idx
  std::vector<int> index_map(nvox, -1);
  for (int m=0; m<M; m++) {
    index_map[masked_indices[m]] = m;
  }

  double h_y = h*h;
  int half_p = (patch_size-1)/2;
  int half_s = search_radius;

  // Precompute offsets for search window
  std::vector<int> dx, dy, dz;
  dx.reserve((2*search_radius+1)*(2*search_radius+1)*(2*search_radius+1));
  dy.reserve(dx.capacity());
  dz.reserve(dx.capacity());

  for (int ddz=-half_s; ddz<=half_s; ddz++){
    for (int ddy=-half_s; ddy<=half_s; ddy++){
      for (int ddx=-half_s; ddx<=half_s; ddx++){
        dx.push_back(ddx); dy.push_back(ddy); dz.push_back(ddz);
      }
    }
  }
  int n_offsets = (int)dx.size();

  // Compute NLM weights for W1 over masked voxels only
  std::vector<arma::uword> Ivec;
  std::vector<arma::uword> Jvec;
  std::vector<double> Vvec;
  Ivec.reserve((size_t)M* (size_t)(n_offsets/5)); // heuristic
  Jvec.reserve(Ivec.capacity());
  Vvec.reserve(Ivec.capacity());

  std::vector<double> center_patch, neigh_patch;
  center_patch.reserve((size_t)patch_size*patch_size*patch_size);
  neigh_patch.reserve(center_patch.size());

  auto idxToCoord = [&](int idx){
    int z_ = idx/(nx*ny);
    int rest = idx%(nx*ny);
    int y_ = rest/nx;
    int x_ = rest%nx;
    return IntegerVector::create(x_, y_, z_);
  };

  arma::vec d_vals(M, arma::fill::zeros);

  for (int m=0; m<M; m++) {
    int i_global = masked_indices[m];
    IntegerVector c = idxToCoord(i_global);
    int x = c[0], y=c[1], z=c[2];
    extract_patch(img_norm, nx, ny, nz, x, y, z, half_p, center_patch);

    for (int o=0; o<n_offsets; o++) {
      int xx = x+dx[o];
      int yy = y+dy[o];
      int zz = z+dz[o];
      if (xx<0 || xx>=nx || yy<0 || yy>=ny || zz<0 || zz>=nz)
        continue;
      int j_global = xx + nx*(yy+ny*zz);
      int j_m = index_map[j_global];
      if (j_m<0) continue; // not in mask

      double dist2=0.0;
      if (j_global==i_global) {
        dist2=0.0;
      } else {
        extract_patch(img_norm, nx, ny, nz, xx, yy, zz, half_p, neigh_patch);
        for (size_t pi=0; pi<center_patch.size(); pi++) {
          double diff = center_patch[pi]-neigh_patch[pi];
          dist2 += diff*diff;
        }
      }
      double kij = std::exp(-dist2/(h_y+1e-15));
      if (kij>1e-8) {
        Ivec.push_back((arma::uword)m);
        Jvec.push_back((arma::uword)j_m);
        Vvec.push_back(kij);
        d_vals[m] += kij;
      }
    }
  }

  {
    arma::umat locs(2, Ivec.size());
    for (size_t i=0; i<Ivec.size(); i++){
      locs(0,i)=Ivec[i]; locs(1,i)=Jvec[i];
    }
    arma::vec vals(Vvec);
    arma::sp_mat K(locs, vals, (arma::uword)M, (arma::uword)M);

    arma::sp_mat W1 = buildW(K, d_vals, use_normalization_free);

    // Compute higher order filters if k>1
    std::vector<arma::sp_mat> W_list; W_list.reserve((size_t)k);
    std::vector<arma::vec> D_list; D_list.reserve((size_t)k);
    W_list.push_back(W1);
    D_list.push_back(d_vals);

    arma::sp_mat prevK = K;
    for (int level=2; level<=k; level++) {
      // K_next = K_l * K_l (hadamard)
      arma::sp_mat K_next = hadamard(prevK, prevK);
      arma::vec d_next(K_next.n_rows, arma::fill::zeros);
      for (auto it = K_next.begin(); it!=K_next.end(); ++it) {
        d_next[it.row()] += (*it);
      }
      arma::sp_mat W_next = buildW(K_next, d_next, use_normalization_free);
      W_list.push_back(W_next);
      D_list.push_back(d_next);
      prevK = K_next;
    }

    // Apply filters to masked voxels
    arma::vec y_vec(M);
    for (int i=0; i<M; i++) {
      y_vec[i] = img_norm[masked_indices[i]];
    }

    arma::vec base_layer = W_list[0]*y_vec;

    std::vector<arma::vec> detail_layers;
    for (int level=1; level<k; level++) {
      arma::vec dl = (W_list[level]*y_vec - W_list[level-1]*y_vec);
      detail_layers.push_back(dl);
    }

    arma::vec high_pass_layer = (y_vec - W_list[k-1]*y_vec);

    // Get mapping params
    auto get_map_params = [&](std::string name, double &a, double &width){
      if (mapping_params.containsElementNamed(name.c_str())) {
        List mp = mapping_params[name];
        a = as<double>(mp["a"]);
        width = as<double>(mp["width"]);
      } else {
        a=10; width=1;
      }
    };

    double a_base, w_base;
    get_map_params("base", a_base, w_base);

    std::vector<double> a_details(k);
    std::vector<double> w_details(k);
    for (int i=0; i<k; i++){
      std::string nm = (i<k-1) ? "detail"+std::to_string(i+1) : "detail"+std::to_string(k);
      double aa, ww;
      get_map_params(nm, aa, ww);
      a_details[i]=aa; w_details[i]=ww;
    }

    for (int i=0; i<M; i++) {
      base_layer[i] = base_map(base_layer[i], a_base, w_base);
    }

    for (int dl=0; dl<(int)detail_layers.size(); dl++) {
      double aa=a_details[dl];
      double ww=w_details[dl];
      for (int i=0; i<M; i++) {
        detail_layers[dl][i] = detail_map(detail_layers[dl][i], aa, ww);
      }
    }

    double a_high=a_details[k-1];
    double w_high=w_details[k-1];
    for (int i=0; i<M; i++){
      high_pass_layer[i] = detail_map(high_pass_layer[i], a_high, w_high);
    }

    // structure mask
    int p = patch_size*patch_size*patch_size;
    arma::vec d_vals_base = D_list[0];
    arma::vec m_vals(M);
    for (int i=0; i<M; i++) {
      double mm = 1.0 - (d_vals_base[i]/(double)p);
      if (mm<0) mm=0;
      if (mm>1) mm=1;
      m_vals[i]=mm;
    }

    // Blend layers
    arma::vec z_vec = base_layer;
    for (int dl=0; dl<(int)detail_layers.size(); dl++) {
      for (int i=0; i<M; i++) {
        z_vec[i] += m_vals[i]*detail_layers[dl][i];
      }
    }
    for (int i=0; i<M; i++) {
      z_vec[i] += m_vals[i]*high_pass_layer[i];
    }

    // Put result back into full volume
    NumericVector out(nvox);
    for (int i=0; i<nvox; i++){
      out[i]=img[i]; // default
    }
    for (int i=0; i<M; i++) {
      out[masked_indices[i]] = z_vec[i]*scale + img_min;
    }
    out.attr("dim") = dims;
    return out;
  }
}
