"""nilearn / scipy reference results and timings.

    NEUROIM2_NL_DIR=/tmp/nl python3 ref_py.py

Writes ref.json plus the arrays compare.py needs. Anything slow is guarded by
NEUROIM2_NL_SCALING=1 (the connected-component scaling probe takes ~15 minutes
on the neuroim2 side, not here).
"""
import json, os, time, warnings
import numpy as np
import nibabel as nib
from scipy import ndimage

warnings.filterwarnings("ignore")
B = os.environ.get("NEUROIM2_NL_DIR", "/tmp/nl")
out = {}


def tm(fn, reps=3):
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter(); fn(); ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


vol = nib.load(f"{B}/vol.nii")
run = nib.load(f"{B}/run.nii")
mask = nib.load(f"{B}/mask.nii")
labels = nib.load(f"{B}/labels.nii")
truth = np.asanyarray(mask.dataobj).astype(bool)

# ---------------------------------------------------------------- smoothing
# The point-spread function is the operator: smooth an impulse and measure the
# FWHM that came out, rather than trusting the argument that went in.
from nilearn.image import (smooth_img, resample_img, resample_to_img, reorder_img,
                           mean_img, index_img, threshold_img, binarize_img,
                           math_img, concat_imgs)

def psf_fwhm(a, spacing, centre):
    prof = a.sum(axis=(1, 2))
    off = (np.arange(len(prof)) - centre) * spacing
    m = prof.sum(); mu = (prof * off).sum() / m
    return float(2 * np.sqrt(2 * np.log(2)) * np.sqrt((prof * (off - mu) ** 2).sum() / m))

imp = np.zeros((41, 41, 41), dtype=np.float32); imp[20, 20, 20] = 1
imp_img = nib.Nifti1Image(imp, np.diag([2., 2, 2, 1]))
out["psf"] = {f"{f}": psf_fwhm(np.asanyarray(smooth_img(imp_img, f).dataobj), 2.0, 20)
              for f in (4.71, 7.06, 9.42)}
out["t_smooth_3d"] = tm(lambda: smooth_img(vol, 6.0))
out["t_smooth_4d"] = tm(lambda: smooth_img(run, 6.0), 1)

# --------------------------------------------------------------- resampling
aff3 = vol.affine.copy(); aff3[:3, :3] = np.diag([-3., 3., 3.])
shape3 = tuple(int(np.ceil(s * 2 / 3)) for s in vol.shape)
for kind, interp in [("linear", "linear"), ("cubic", "continuous")]:
    img = resample_img(vol, aff3, shape3, interpolation=interp,
                       force_resample=True, copy_header=True)
    nib.save(img, f"{B}/py_down_{kind}.nii")
out["t_resample"] = tm(lambda: resample_img(vol, aff3, shape3, interpolation="linear",
                                            force_resample=True, copy_header=True))
same = resample_img(vol, vol.affine, vol.shape, interpolation="continuous",
                    force_resample=True, copy_header=True)
out["identity_resample_maxdiff"] = float(
    np.abs(np.asanyarray(vol.dataobj) - np.asanyarray(same.dataobj)).max())

dlab = resample_img(labels, aff3, shape3, interpolation="nearest",
                    force_resample=True, copy_header=True)
dl = np.asanyarray(dlab.dataobj)
out["nearest_labels_in"] = int(len(np.unique(np.asanyarray(labels.dataobj))))
out["nearest_labels_out"] = int(len(np.unique(dl)))

# ------------------------------------------------------------- reorientation
odd = nib.load(f"{B}/odd.nii")
ro = reorder_img(odd)
a0 = np.asanyarray(odd.dataobj); a1 = np.asanyarray(ro.dataobj)
out["reorder"] = {
    "in_shape": list(odd.shape), "out_shape": list(ro.shape),
    "in_axcodes": "".join(nib.aff2axcodes(odd.affine)),
    "out_axcodes": "".join(nib.aff2axcodes(ro.affine)),
    "nonzero_in": int((a0 != 0).sum()), "nonzero_out": int((a1 != 0).sum()),
    "values_are_a_permutation": bool(np.array_equal(np.sort(a0.ravel()), np.sort(a1.ravel()))),
    "seconds": tm(lambda: reorder_img(odd)),
}

# ------------------------------------------------- connected components
a = np.asanyarray(vol.dataobj)
binv = a > 110
STRUCT = {
    "6-connect": ndimage.generate_binary_structure(3, 1),
    "18-connect": np.array([[[0, 1, 0], [1, 1, 1], [0, 1, 0]],
                            [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
                            [[0, 1, 0], [1, 1, 1], [0, 1, 0]]]),
    "26-connect": ndimage.generate_binary_structure(3, 3),
}
out["cc"] = {}
for name, st in STRUCT.items():
    lab, n = ndimage.label(binv, structure=st)
    out["cc"][name] = {"n": int(n), "largest": int(np.bincount(lab.ravel())[1:].max()),
                       "seconds": tm(lambda st=st: ndimage.label(binv, structure=st))}
out["cc_in_mask"] = int(binv.sum())

if os.environ.get("NEUROIM2_NL_SCALING"):
    out["cc_scaling"] = {}
    st = STRUCT["26-connect"]
    for name in ("small", "typical", "hires"):
        d = np.asanyarray(nib.load(f"{B}/scale_{name}.nii").dataobj) > 0.9
        lab, n = ndimage.label(d, structure=st)
        out["cc_scaling"][name] = {"in_mask": int(d.sum()), "n": int(n),
                                   "seconds": tm(lambda d=d, st=st: ndimage.label(d, structure=st))}

# ------------------------------------------------------------------ masking
from nilearn.masking import compute_epi_mask


def dice(x):
    return float(2 * (x & truth).sum() / (x.sum() + truth.sum()))


em = compute_epi_mask(run)
ea = np.asanyarray(em.dataobj).astype(bool)
out["epi_mask"] = {"voxels": int(ea.sum()), "dice": dice(ea),
                   "seconds": tm(lambda: compute_epi_mask(run))}

# ------------------------------------------------- spheres, parcels, masking
from nilearn.maskers import NiftiSpheresMasker, NiftiLabelsMasker, NiftiMasker

world = np.loadtxt(f"{B}/sphere_centres_world.txt")
out["spheres"] = {}
for rad in (6.0, 8.0):
    m = NiftiSpheresMasker([tuple(w) for w in world], radius=rad, mask_img=mask,
                           allow_overlap=True, standardize=False)
    sig = m.fit_transform(run)
    np.save(f"{B}/py_spheres_{int(rad)}.npy", sig)
    out["spheres"][str(int(rad))] = {"shape": list(sig.shape),
                                     "seconds": tm(lambda m=m: m.fit_transform(run), 1)}

lm = NiftiLabelsMasker(labels_img=labels, standardize=False)
ls = lm.fit_transform(run)
np.save(f"{B}/py_labels.npy", ls)
out["labels_signal"] = {"shape": list(ls.shape), "seconds": tm(lambda: lm.fit_transform(run), 1)}

nmk = NiftiMasker(mask_img=mask, standardize=False)
nmk.fit_transform(run)
out["masked_extract"] = {"seconds": tm(lambda: nmk.fit_transform(run), 1)}

# ------------------------------------------------- searchlight neighbourhoods
from sklearn import neighbors
idx = np.argwhere(truth)
wc = nib.affines.apply_affine(mask.affine, idx)
out["searchlight"] = {}
for rad in (6.0, 8.0):
    def build(rad=rad):
        clf = neighbors.NearestNeighbors(radius=rad).fit(wc)
        return clf.radius_neighbors_graph(wc, mode="connectivity")
    A = build()
    sizes = np.diff(A.indptr)
    out["searchlight"][str(int(rad))] = {
        "centres": int(len(idx)), "mean": float(sizes.mean()),
        "median": int(np.median(sizes)), "seconds": tm(build, 1)}

# ----------------------------------------------------------------- utilities
np.save(f"{B}/py_mean.npy", np.asanyarray(mean_img(run, copy_header=True).dataobj))
sub = index_img(run, slice(0, 20))
out["util"] = {
    "mean_img": tm(lambda: mean_img(run, copy_header=True)),
    "index_img": tm(lambda: index_img(run, 30)),
    "threshold_img": tm(lambda: threshold_img(vol, 110, copy_header=True,
                                              cluster_threshold=0, two_sided=False)),
    "binarize_img": tm(lambda: binarize_img(vol, 110, copy_header=True)),
    "math_img": tm(lambda: math_img("img*2+1", img=vol)),
    "concat_imgs": tm(lambda: concat_imgs([sub, sub, sub]), 1),
}

with open(f"{B}/ref.json", "w") as f:
    json.dump(out, f, indent=1)
print(f"wrote {B}/ref.json")
