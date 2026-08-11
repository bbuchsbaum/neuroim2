"""Generate a battery of NIfTI conformance probes with nibabel (the reference)."""
import os, struct, gzip
import numpy as np, nibabel as nib

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "probes")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(7)
AFF = np.array([[-2., 0, 0, 90.], [0, 2., 0, -126.], [0, 0, 2., -72.], [0, 0, 0, 1.]])
SHAPE = (8, 9, 7)

manifest = {}

def save(name, img):
    p = f"{OUT}/{name}.nii"
    img.to_filename(p)
    manifest[name] = p
    return p

# --- data types ---
for dt, code in [(np.uint8, 'uint8'), (np.int16, 'int16'), (np.int32, 'int32'),
                 (np.float32, 'float32'), (np.float64, 'float64'),
                 (np.uint16, 'uint16'), (np.uint32, 'uint32'),
                 (np.int64, 'int64'), (np.uint64, 'uint64'),
                 (np.int8, 'int8')]:
    info = np.iinfo(dt) if np.issubdtype(dt, np.integer) else None
    if info is not None:
        lo = max(info.min, -10**6); hi = min(info.max, 10**6)
        d = rng.integers(lo, hi, size=SHAPE).astype(dt)
        # pin extremes so signedness errors are visible
        d.flat[0] = info.min; d.flat[1] = info.max
    else:
        d = rng.normal(size=SHAPE).astype(dt)
    img = nib.Nifti1Image(d, AFF, dtype=dt)
    try:
        save(f"dtype_{code}", img)
    except Exception as e:
        print("skip", code, e)

# complex + RGB
d = (rng.normal(size=SHAPE) + 1j*rng.normal(size=SHAPE)).astype(np.complex64)
save("dtype_complex64", nib.Nifti1Image(d, AFF, dtype=np.complex64))
rgb = np.dtype([('R', 'u1'), ('G', 'u1'), ('B', 'u1')])
d = np.zeros(SHAPE, dtype=rgb)
d['R'] = rng.integers(0, 255, SHAPE); d['G'] = 7; d['B'] = 200
save("dtype_rgb24", nib.Nifti1Image(d, AFF, dtype=rgb))

# --- scaling conventions ---
base = rng.integers(-1000, 1000, size=SHAPE).astype(np.int16)
for name, slope, inter in [("scl_2_5", 2.0, 5.0), ("scl_zero", 0.0, 0.0),
                           ("scl_nan", np.nan, np.nan), ("scl_neg", -1.5, 3.0),
                           ("scl_inter_only", 1.0, -100.0)]:
    img = nib.Nifti1Image(base, AFF, dtype=np.int16)
    p = f"{OUT}/{name}.nii"
    img.to_filename(p)
    # patch the raw header so nibabel's writer cannot normalise our values away
    with open(p, 'r+b') as f:
        f.seek(112); f.write(struct.pack('<ff', slope, inter))
    manifest[name] = p

# --- affine provenance ---
img = nib.Nifti1Image(base, AFF, dtype=np.int16)
p = f"{OUT}/sform_only.nii"; img.to_filename(p)
with open(p, 'r+b') as f:
    f.seek(252); f.write(struct.pack('<hh', 0, 4))       # qform_code=0, sform_code=4
manifest["sform_only"] = p

p = f"{OUT}/qform_only.nii"; img.to_filename(p)
with open(p, 'r+b') as f:
    f.seek(252); f.write(struct.pack('<hh', 1, 0))       # qform_code=1, sform_code=0
    f.seek(280); f.write(struct.pack('<12f', *([0.0]*12)))  # zero the srow rows
manifest["qform_only"] = p

p = f"{OUT}/no_form.nii"; img.to_filename(p)
with open(p, 'r+b') as f:
    f.seek(252); f.write(struct.pack('<hh', 0, 0))
    f.seek(280); f.write(struct.pack('<12f', *([0.0]*12)))
manifest["no_form"] = p

# qform with qfac = -1 (left-handed)
aff_lh = AFF.copy(); aff_lh[:3, :3] = np.array([[2., 0, 0], [0, 2., 0], [0, 0, 2.]])
img = nib.Nifti1Image(base, aff_lh, dtype=np.int16)
img.header.set_qform(aff_lh, code=1); img.header.set_sform(np.zeros((4, 4)), code=0)
p = f"{OUT}/qfac_neg.nii"
img.to_filename(p)
with open(p, 'r+b') as f:
    f.seek(76); f.write(struct.pack('<f', -1.0))          # pixdim[0] = -1
    f.seek(252); f.write(struct.pack('<hh', 1, 0))
    f.seek(280); f.write(struct.pack('<12f', *([0.0]*12)))
manifest["qfac_neg"] = p

# --- endianness ---
be_hdr = nib.Nifti1Header(endianness='>')
be_hdr.set_data_shape(SHAPE); be_hdr.set_data_dtype(np.int16)
be_hdr.set_sform(AFF, code=1); be_hdr.set_qform(AFF, code=1)
p = f"{OUT}/bigendian.nii"
nib.Nifti1Image(base, AFF, header=be_hdr).to_filename(p)
manifest["bigendian"] = p

# --- two-file .hdr/.img pair ---
pair = nib.Nifti1Pair(base, AFF)
pair.to_filename(f"{OUT}/pair.hdr")
manifest["pair"] = f"{OUT}/pair.hdr"

# --- ANALYZE 7.5 ---
ana = nib.AnalyzeImage(base, AFF)
ana.to_filename(f"{OUT}/analyze.hdr")
manifest["analyze"] = f"{OUT}/analyze.hdr"

# --- NIfTI-2 ---
n2 = nib.Nifti2Image(base, AFF)
n2.to_filename(f"{OUT}/nifti2.nii")
manifest["nifti2"] = f"{OUT}/nifti2.nii"

# --- gzip ---
nib.Nifti1Image(base, AFF, dtype=np.int16).to_filename(f"{OUT}/gz.nii.gz")
manifest["gz"] = f"{OUT}/gz.nii.gz"

# --- non-finite float payload ---
d = rng.normal(size=SHAPE).astype(np.float32)
d.flat[0] = np.nan; d.flat[1] = np.inf; d.flat[2] = -np.inf
nib.Nifti1Image(d, AFF, dtype=np.float32).to_filename(f"{OUT}/nonfinite.nii")
manifest["nonfinite"] = f"{OUT}/nonfinite.nii"

# --- 5D / 4D-with-singleton ---
d5 = rng.integers(0, 100, size=(6, 7, 5, 1, 3)).astype(np.int16)
nib.Nifti1Image(d5, AFF, dtype=np.int16).to_filename(f"{OUT}/five_d.nii")
manifest["five_d"] = f"{OUT}/five_d.nii"

d4s = rng.integers(0, 100, size=(6, 7, 5, 1)).astype(np.int16)
nib.Nifti1Image(d4s, AFF, dtype=np.int16).to_filename(f"{OUT}/four_d_singleton.nii")
manifest["four_d_singleton"] = f"{OUT}/four_d_singleton.nii"

# --- header extension present (vox_offset > 352) ---
img = nib.Nifti1Image(base, AFF, dtype=np.int16)
img.header.extensions.append(nib.nifti1.Nifti1Extension(6, b'a comment extension'))
img.to_filename(f"{OUT}/with_extension.nii")
manifest["with_extension"] = f"{OUT}/with_extension.nii"

# --- expected values, computed by nibabel ---
import json
exp = {}
for name, path in manifest.items():
    try:
        im = nib.load(path)
        arr = np.asanyarray(im.dataobj)
        if arr.dtype.names:   # RGB struct
            summ = {"kind": "struct", "first": [int(x) for x in arr.flat[0]]}
        elif np.iscomplexobj(arr):
            summ = {"kind": "complex", "first_re": float(arr.flat[0].real),
                    "first_im": float(arr.flat[0].imag)}
        else:
            a = arr.astype(np.float64)
            fin = a[np.isfinite(a)]
            summ = {"kind": "numeric",
                    "shape": list(arr.shape),
                    "first": float(a.flat[0]), "second": float(a.flat[1]),
                    "min": float(fin.min()) if fin.size else None,
                    "max": float(fin.max()) if fin.size else None,
                    "sum": float(fin.sum()) if fin.size else None,
                    "n_nan": int(np.isnan(a).sum()), "n_inf": int(np.isinf(a).sum())}
        exp[name] = {"path": path, "affine": im.affine.tolist(), "data": summ}
    except Exception as e:
        exp[name] = {"path": path, "error": f"{type(e).__name__}: {e}"}

with open(f"{OUT}/expected.json", "w") as f:
    json.dump(exp, f, indent=1)
print("probes:", len(manifest))
for k in sorted(exp):
    if "error" in exp[k]:
        print("  nibabel could not load", k, exp[k]["error"])
