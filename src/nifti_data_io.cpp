// Compiled binary payload I/O for NIfTI (and other raw-block) images.
//
// R's readBin/writeBin go through the connection layer and convert one element
// at a time; measured on 7.2M int16 values that is roughly 9x the cost of the
// equivalent fread-and-convert loop. Everything here is a straight
// file -> double (or double -> file) pass with the element conversion inlined,
// which is what the read and write paths in R/binary_io.R now call.
//
// Reading always produces double. That is not only for speed: R's integer type
// reserves INT_MIN for NA, so an int32 voxel legitimately holding -2147483648
// becomes NA if it is routed through an R integer vector.

#include <Rcpp.h>
#include <zlib.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

// NIfTI datatype codes we handle as a raw block of fixed-width elements.
static const int NI_UINT8 = 2;
static const int NI_INT16 = 4;
static const int NI_INT32 = 8;
static const int NI_FLOAT32 = 16;
static const int NI_FLOAT64 = 64;
static const int NI_INT8 = 256;
static const int NI_UINT16 = 512;
static const int NI_UINT32 = 768;
static const int NI_INT64 = 1024;
static const int NI_UINT64 = 1280;

static int dtype_width(int code) {
  switch (code) {
    case NI_UINT8:
    case NI_INT8:
      return 1;
    case NI_INT16:
    case NI_UINT16:
      return 2;
    case NI_INT32:
    case NI_UINT32:
    case NI_FLOAT32:
      return 4;
    case NI_FLOAT64:
    case NI_INT64:
    case NI_UINT64:
      return 8;
    default:
      return 0;
  }
}

static void swap_in_place(unsigned char* p, size_t n_elem, int width) {
  if (width < 2) return;
  for (size_t i = 0; i < n_elem; ++i) {
    unsigned char* e = p + (size_t) i * width;
    for (int a = 0, b = width - 1; a < b; ++a, --b) std::swap(e[a], e[b]);
  }
}

// Decode `n` elements of `code` from `src` into `dst`.
static void decode_block(const unsigned char* src, double* dst, size_t n, int code) {
  switch (code) {
    case NI_UINT8: {
      for (size_t i = 0; i < n; ++i) dst[i] = (double) src[i];
      break;
    }
    case NI_INT8: {
      const int8_t* p = reinterpret_cast<const int8_t*>(src);
      for (size_t i = 0; i < n; ++i) dst[i] = (double) p[i];
      break;
    }
    case NI_INT16: {
      int16_t v;
      for (size_t i = 0; i < n; ++i) { std::memcpy(&v, src + 2 * i, 2); dst[i] = (double) v; }
      break;
    }
    case NI_UINT16: {
      uint16_t v;
      for (size_t i = 0; i < n; ++i) { std::memcpy(&v, src + 2 * i, 2); dst[i] = (double) v; }
      break;
    }
    case NI_INT32: {
      int32_t v;
      for (size_t i = 0; i < n; ++i) { std::memcpy(&v, src + 4 * i, 4); dst[i] = (double) v; }
      break;
    }
    case NI_UINT32: {
      uint32_t v;
      for (size_t i = 0; i < n; ++i) { std::memcpy(&v, src + 4 * i, 4); dst[i] = (double) v; }
      break;
    }
    case NI_FLOAT32: {
      float v;
      for (size_t i = 0; i < n; ++i) { std::memcpy(&v, src + 4 * i, 4); dst[i] = (double) v; }
      break;
    }
    case NI_FLOAT64: {
      double v;
      for (size_t i = 0; i < n; ++i) { std::memcpy(&v, src + 8 * i, 8); dst[i] = v; }
      break;
    }
    case NI_INT64: {
      int64_t v;
      for (size_t i = 0; i < n; ++i) { std::memcpy(&v, src + 8 * i, 8); dst[i] = (double) v; }
      break;
    }
    case NI_UINT64: {
      uint64_t v;
      for (size_t i = 0; i < n; ++i) { std::memcpy(&v, src + 8 * i, 8); dst[i] = (double) v; }
      break;
    }
    default:
      stop("decode_block: unhandled datatype code %d", code);
  }
}

// Encode `n` doubles into `dst` as `code`. Integer targets round half-to-even
// (matching R's round() and numpy's rint) and saturate at the type bounds;
// a non-finite value has no integer representation and is written as 0.
static void encode_block(const double* src, unsigned char* dst, size_t n, int code) {
  switch (code) {
    case NI_FLOAT64: {
      std::memcpy(dst, src, n * sizeof(double));
      break;
    }
    case NI_FLOAT32: {
      for (size_t i = 0; i < n; ++i) {
        float v = (float) src[i];
        std::memcpy(dst + 4 * i, &v, 4);
      }
      break;
    }
    default:
      break;
  }
  if (code == NI_FLOAT32 || code == NI_FLOAT64) return;

  double lo, hi;
  switch (code) {
    case NI_UINT8:  lo = 0;                     hi = 255;                  break;
    case NI_INT8:   lo = -128;                  hi = 127;                  break;
    case NI_INT16:  lo = -32768;                hi = 32767;                break;
    case NI_UINT16: lo = 0;                     hi = 65535;                break;
    case NI_INT32:  lo = -2147483648.0;         hi = 2147483647.0;         break;
    case NI_UINT32: lo = 0;                     hi = 4294967295.0;         break;
    // The exact bounds of int64/uint64 are not representable in double; use
    // the largest magnitudes that are, so the saturation itself is exact.
    case NI_INT64:  lo = -9223372036854775808.0; hi = 9223372036854774784.0; break;
    case NI_UINT64: lo = 0;                      hi = 18446744073709549568.0; break;
    default: stop("encode_block: unhandled datatype code %d", code);
  }

  for (size_t i = 0; i < n; ++i) {
    double x = src[i];
    double r;
    if (!R_FINITE(x)) {
      r = 0.0;
    } else {
      r = std::nearbyint(x);
      if (r < lo) r = lo;
      if (r > hi) r = hi;
    }
    switch (code) {
      case NI_UINT8:  { uint8_t  v = (uint8_t)  r; dst[i] = v; break; }
      case NI_INT8:   { int8_t   v = (int8_t)   r; std::memcpy(dst + i, &v, 1); break; }
      case NI_INT16:  { int16_t  v = (int16_t)  r; std::memcpy(dst + 2 * i, &v, 2); break; }
      case NI_UINT16: { uint16_t v = (uint16_t) r; std::memcpy(dst + 2 * i, &v, 2); break; }
      case NI_INT32:  { int32_t  v = (int32_t)  r; std::memcpy(dst + 4 * i, &v, 4); break; }
      case NI_UINT32: { uint32_t v = (uint32_t) r; std::memcpy(dst + 4 * i, &v, 4); break; }
      case NI_INT64:  { int64_t  v = (int64_t)  r; std::memcpy(dst + 8 * i, &v, 8); break; }
      case NI_UINT64: { uint64_t v = (uint64_t) r; std::memcpy(dst + 8 * i, &v, 8); break; }
    }
  }
}

// A file handle that reads either a plain or a gzipped stream.
namespace {
struct InStream {
  FILE* raw = nullptr;
  gzFile gz = nullptr;

  bool open(const std::string& path, bool gzipped) {
    if (gzipped) {
      gz = gzopen(path.c_str(), "rb");
      return gz != nullptr;
    }
    raw = fopen(path.c_str(), "rb");
    return raw != nullptr;
  }
  // Advance to a byte offset. gzseek decompresses up to that point, which is
  // unavoidable for a gzip stream and is what the R path did too.
  bool seek(double offset) {
    if (gz) return gzseek(gz, (z_off_t) offset, SEEK_SET) >= 0;
    return fseeko(raw, (off_t) offset, SEEK_SET) == 0;
  }
  size_t read(void* buf, size_t nbytes) {
    if (gz) {
      int got = gzread(gz, buf, (unsigned int) nbytes);
      return got < 0 ? 0 : (size_t) got;
    }
    return fread(buf, 1, nbytes, raw);
  }
  void close() {
    if (gz) { gzclose(gz); gz = nullptr; }
    if (raw) { fclose(raw); raw = nullptr; }
  }
  ~InStream() { close(); }
};
}  // namespace

//' Read a raw block of image data into a double vector
//'
//' @param path file to read from
//' @param offset byte offset of the first element
//' @param n number of elements to read
//' @param dtype_code NIfTI datatype code of the stored elements
//' @param swap TRUE when the file's byte order differs from the platform's
//' @param gzipped TRUE when the file is a gzip stream
//' @return a numeric vector of length \code{n}
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector nifti_read_data_cpp(std::string path, double offset, double n,
                                  int dtype_code, bool swap, bool gzipped) {
  if (!R_FINITE(n) || n < 0) stop("'n' must be a non-negative finite count");
  if (!R_FINITE(offset) || offset < 0) stop("'offset' must be a non-negative finite byte offset");
  int width = dtype_width(dtype_code);
  if (width == 0) stop("nifti_read_data_cpp: unhandled datatype code %d", dtype_code);
  if (n > (double) R_XLEN_T_MAX) stop("requested element count exceeds R's vector limit");

  R_xlen_t N = (R_xlen_t) n;
  NumericVector out(N);
  if (N == 0) return out;

  InStream in;
  if (!in.open(path, gzipped)) stop("cannot open '%s' for reading", path.c_str());
  if (!in.seek(offset))
    stop("cannot seek to byte %.0f in '%s'; the file is shorter than its header claims",
         offset, path.c_str());

  const size_t CHUNK = 1u << 20;  // elements
  std::vector<unsigned char> buf(CHUNK * (size_t) width);
  double* o = REAL(out);
  R_xlen_t done = 0;

  while (done < N) {
    size_t want = (size_t) std::min<R_xlen_t>((R_xlen_t) CHUNK, N - done);
    size_t want_bytes = want * (size_t) width;
    size_t got_bytes = 0;
    // A short read is normal for a gzip stream near an internal buffer edge,
    // so keep asking until the request is satisfied or the file really ends.
    while (got_bytes < want_bytes) {
      size_t r = in.read(buf.data() + got_bytes, want_bytes - got_bytes);
      if (r == 0) break;
      got_bytes += r;
    }
    if (got_bytes < want_bytes) {
      double have = (double) done + (double) (got_bytes / (size_t) width);
      in.close();
      stop("'%s' is truncated: the header describes %.0f elements of %d bytes "
           "beginning at byte %.0f (%.0f bytes of data), but the file supplies only %.0f",
           path.c_str(), n, width, offset, n * width, have);
    }
    if (swap) swap_in_place(buf.data(), want, width);
    decode_block(buf.data(), o + done, want, dtype_code);
    done += (R_xlen_t) want;
  }
  in.close();
  return out;
}

//' Write a header block followed by image data
//'
//' @param path file to create
//' @param header raw vector written verbatim before the data
//' @param data numeric values to encode
//' @param dtype_code NIfTI datatype code to encode as
//' @param slope,inter scaling to invert before encoding, i.e. the stored value
//'   is \code{(x - inter) / slope}
//' @param swap TRUE to write in the opposite byte order from the platform's
//' @param gzipped TRUE to gzip the output
//' @return the number of data elements written, invisibly
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double nifti_write_data_cpp(std::string path, RawVector header, NumericVector data,
                            int dtype_code, double slope, double inter,
                            bool swap, bool gzipped) {
  int width = dtype_width(dtype_code);
  if (width == 0) stop("nifti_write_data_cpp: unhandled datatype code %d", dtype_code);
  if (!R_FINITE(slope) || slope == 0.0) slope = 1.0;
  if (!R_FINITE(inter)) inter = 0.0;

  FILE* raw = nullptr;
  gzFile gz = nullptr;
  if (gzipped) {
    gz = gzopen(path.c_str(), "wb");
    if (!gz) stop("cannot open '%s' for writing", path.c_str());
  } else {
    raw = fopen(path.c_str(), "wb");
    if (!raw) stop("cannot open '%s' for writing", path.c_str());
  }

  struct Closer {
    FILE* raw; gzFile gz;
    ~Closer() { if (gz) gzclose(gz); if (raw) fclose(raw); }
  } closer{raw, gz};

  auto put = [&](const void* p, size_t nbytes) -> bool {
    if (nbytes == 0) return true;
    if (gz) return gzwrite(gz, p, (unsigned int) nbytes) == (int) nbytes;
    return fwrite(p, 1, nbytes, raw) == nbytes;
  };

  if (header.size() > 0) {
    if (!put(&header[0], (size_t) header.size()))
      stop("failed writing the header block of '%s'", path.c_str());
  }

  R_xlen_t N = data.size();
  const double* d = REAL(data);
  const bool rescale = (slope != 1.0 || inter != 0.0);
  const size_t CHUNK = 1u << 20;
  std::vector<unsigned char> buf(CHUNK * (size_t) width);
  std::vector<double> staged;
  if (rescale) staged.resize(CHUNK);

  R_xlen_t done = 0;
  while (done < N) {
    size_t take = (size_t) std::min<R_xlen_t>((R_xlen_t) CHUNK, N - done);
    const double* src = d + done;
    if (rescale) {
      for (size_t i = 0; i < take; ++i) staged[i] = (src[i] - inter) / slope;
      src = staged.data();
    }
    encode_block(src, buf.data(), take, dtype_code);
    if (swap) swap_in_place(buf.data(), take, width);
    if (!put(buf.data(), take * (size_t) width))
      stop("failed writing image data to '%s'; the device may be full", path.c_str());
    done += (R_xlen_t) take;
  }
  return (double) N;
}

//' Gather a set of volumes from a 4-D image file
//'
//' Reads whole volumes rather than an arbitrary index set, so a contiguous
//' request costs one sequential pass. Returns a
//' \code{prod(dim[1:3]) x length(vols)} matrix, which is the layout
//' \code{DenseNeuroVec} wants, so no transpose is needed on either side.
//'
//' @param path file to read from
//' @param offset byte offset of the first element of the first volume
//' @param nels voxels per volume
//' @param vols 1-based volume indices, in the order they should appear
//' @param dtype_code NIfTI datatype code of the stored elements
//' @param swap TRUE when the file's byte order differs from the platform's
//' @param gzipped TRUE when the file is a gzip stream
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix nifti_read_volumes_cpp(std::string path, double offset, double nels,
                                     NumericVector vols, int dtype_code,
                                     bool swap, bool gzipped) {
  int width = dtype_width(dtype_code);
  if (width == 0) stop("nifti_read_volumes_cpp: unhandled datatype code %d", dtype_code);
  if (!R_FINITE(nels) || nels < 0) stop("'nels' must be a non-negative finite count");
  R_xlen_t V = vols.size();
  R_xlen_t NE = (R_xlen_t) nels;
  if (V > 0 && nels > 0 && (double) R_XLEN_T_MAX / nels < (double) V)
    stop("requested volumes exceed R's vector limit");

  NumericMatrix out(NE, V);
  if (V == 0 || NE == 0) return out;

  InStream in;
  if (!in.open(path, gzipped)) stop("cannot open '%s' for reading", path.c_str());

  const size_t CHUNK = 1u << 20;
  std::vector<unsigned char> buf(CHUNK * (size_t) width);
  double* o = REAL(out);
  double expect_pos = -1.0;  // byte offset the stream is already at, if known

  for (R_xlen_t v = 0; v < V; ++v) {
    double idx = vols[v];
    if (!R_FINITE(idx) || idx < 1) stop("volume indices must be finite and >= 1");
    double vol_off = offset + (idx - 1.0) * nels * width;
    if (vol_off != expect_pos) {
      if (!in.seek(vol_off))
        stop("cannot seek to byte %.0f in '%s'; the file is shorter than its header claims",
             vol_off, path.c_str());
    }
    R_xlen_t done = 0;
    double* dst = o + (R_xlen_t) v * NE;
    while (done < NE) {
      size_t want = (size_t) std::min<R_xlen_t>((R_xlen_t) CHUNK, NE - done);
      size_t want_bytes = want * (size_t) width;
      size_t got_bytes = 0;
      while (got_bytes < want_bytes) {
        size_t r = in.read(buf.data() + got_bytes, want_bytes - got_bytes);
        if (r == 0) break;
        got_bytes += r;
      }
      if (got_bytes < want_bytes) {
        in.close();
        stop("'%s' is truncated: volume %.0f needs %.0f bytes beginning at byte %.0f, "
             "but the file supplies fewer", path.c_str(), idx, nels * width, vol_off);
      }
      if (swap) swap_in_place(buf.data(), want, width);
      decode_block(buf.data(), dst + done, want, dtype_code);
      done += (R_xlen_t) want;
    }
    expect_pos = vol_off + nels * width;
  }
  in.close();
  return out;
}
