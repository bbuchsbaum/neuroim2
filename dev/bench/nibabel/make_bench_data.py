"""Generate the timing inputs used by bench_r.R and bench_py.py.

    NEUROIM2_BENCH_DIR=/tmp/bench python3 make_bench_data.py
"""
import os
import numpy as np
import nibabel as nib

B = os.environ.get("NEUROIM2_BENCH_DIR", "/tmp/bench")
os.makedirs(B, exist_ok=True)
rng = np.random.default_rng(0)

# 3D, roughly a 1 mm structural
aff = np.array([[-2., 0, 0, 90.], [0, 2., 0, -126.], [0, 0, 2., -72.], [0, 0, 0, 1.]])
d3 = rng.normal(size=(182, 218, 182)).astype(np.float32)
nib.save(nib.Nifti1Image(d3, aff), f"{B}/vol3d_f32.nii")
nib.save(nib.Nifti1Image(d3, aff), f"{B}/vol3d_f32.nii.gz")

d3i = (rng.normal(size=(182, 218, 182)) * 100).astype(np.int16)
nib.save(nib.Nifti1Image(d3i, aff, dtype=np.int16), f"{B}/vol3d_i16.nii")

# 4D, an ordinary single-run fMRI acquisition
aff4 = np.array([[-3., 0, 0, 96.], [0, 3., 0, -96.], [0, 0, 3.3, -54.], [0, 0, 0, 1.]])
d4 = (rng.normal(size=(64, 64, 36, 200)) * 100).astype(np.int16)
img4 = nib.Nifti1Image(d4, aff4, dtype=np.int16)
img4.header["pixdim"][4] = 2.0            # TR
img4.header.set_xyzt_units("mm", "sec")
nib.save(img4, f"{B}/fmri_i16.nii")
nib.save(img4, f"{B}/fmri_i16.nii.gz")

for f in sorted(os.listdir(B)):
    if f.endswith((".nii", ".nii.gz")):
        print(f"{f:<24}{os.path.getsize(os.path.join(B, f)) // 1024:>10} KB")
