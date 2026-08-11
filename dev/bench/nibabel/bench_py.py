import time, statistics, tempfile, os
import numpy as np, nibabel as nib

B = os.environ.get("NEUROIM2_BENCH_DIR", "/tmp/bench")

def tm(label, fn, reps=3):
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter(); fn(); ts.append(time.perf_counter() - t0)
    print(f"{label:<46} {statistics.median(ts):8.3f} s")
    return statistics.median(ts)

tm("header  3D .nii (nib.load)",          lambda: nib.load(f"{B}/vol3d_i16.nii").header)
tm("header  4D .nii (nib.load)",          lambda: nib.load(f"{B}/fmri_i16.nii").header)
tm("read    3D int16 .nii (get_fdata)",   lambda: nib.load(f"{B}/vol3d_i16.nii").get_fdata())
tm("read    3D float32 .nii (get_fdata)", lambda: nib.load(f"{B}/vol3d_f32.nii").get_fdata())
tm("read    3D f32 .nii.gz (get_fdata)",  lambda: nib.load(f"{B}/vol3d_f32.nii.gz").get_fdata())
tm("read    4D int16 .nii (get_fdata)",   lambda: nib.load(f"{B}/fmri_i16.nii").get_fdata(), reps=3)
tm("read    4D int16 .nii.gz (get_fdata)",lambda: nib.load(f"{B}/fmri_i16.nii.gz").get_fdata(), reps=2)
tm("read    4D subvolume k=100 (dataobj)",lambda: np.asarray(nib.load(f"{B}/fmri_i16.nii").dataobj[..., 99]))

print()
tm("read    3D int16 native (dataobj)",   lambda: np.asanyarray(nib.load(f"{B}/vol3d_i16.nii").dataobj))
tm("read    4D int16 native (dataobj)",   lambda: np.asanyarray(nib.load(f"{B}/fmri_i16.nii").dataobj))

img = nib.load(f"{B}/vol3d_f32.nii")
d = img.get_fdata()
out = nib.Nifti1Image(d.astype(np.float32), img.affine)
tm("write   3D float32 .nii",             lambda: nib.save(out, tempfile.mktemp(suffix=".nii")))
tm("write   3D float32 .nii.gz",          lambda: nib.save(out, tempfile.mktemp(suffix=".nii.gz")), reps=2)

print()
arr = np.asanyarray(nib.load(f"{B}/fmri_i16.nii").dataobj)
flat = arr.reshape(-1, arr.shape[-1])
idx = np.linspace(0, flat.shape[0]-1, 1000).astype(int)
tm("series  1000 voxels x 200 tp (in-mem)", lambda: flat[idx, :].astype(np.float64), reps=5)

print("\n-- lazy access (ArrayProxy, no full load) --")
tm("lazy    open 4D (nib.load)",           lambda: nib.load(f"{B}/fmri_i16.nii"))
def lazy_series():
    im = nib.load(f"{B}/fmri_i16.nii")
    sh = im.shape
    ijk = np.unravel_index(idx, sh[:3])
    # nibabel fancy-index one voxel at a time is the honest comparison for scattered voxels
    return np.stack([im.dataobj[i, j, k, :] for i, j, k in zip(*ijk)])
tm("lazy    series 1000 voxels x 200 tp", lazy_series, reps=3)
