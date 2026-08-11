"""Shared inputs for the nilearn comparison.

    NEUROIM2_NL_DIR=/tmp/nl python3 make_data.py

Everything both sides read comes from these files, so neither library is
handed a different array than the other.
"""
import os
import numpy as np
import nibabel as nib

B = os.environ.get("NEUROIM2_NL_DIR", "/tmp/nl")
os.makedirs(B, exist_ok=True)
rng = np.random.default_rng(11)

AFF = np.array([[-2., 0, 0, 60.], [0, 2., 0, -72.], [0, 0, 2., -56.], [0, 0, 0, 1.]])
SHAPE = (60, 72, 56)

# An ellipsoidal "brain" with signal inside and near-zero outside, so a mask
# estimated from the data has a ground truth to be scored against.
zz, yy, xx = np.meshgrid(np.arange(56), np.arange(72), np.arange(60), indexing="ij")
r = np.sqrt(((xx - 30) / 26) ** 2 + ((yy - 36) / 32) ** 2 + ((zz - 28) / 24) ** 2)
brain = (r < 1.0).T

d = (rng.normal(size=SHAPE) * 20 + 100).astype(np.float32)
d[~brain] = rng.normal(size=(~brain).sum()) * 3
nib.save(nib.Nifti1Image(d, AFF), f"{B}/vol.nii")
nib.save(nib.Nifti1Image(brain.astype(np.uint8), AFF, dtype=np.uint8), f"{B}/mask.nii")

T = 60
d4 = (rng.normal(size=SHAPE + (T,)) * 10 + 100).astype(np.float32)
d4[~brain] = 0
img4 = nib.Nifti1Image(d4, AFF)
img4.header["pixdim"][4] = 2.0
img4.header.set_xyzt_units("mm", "sec")
nib.save(img4, f"{B}/run.nii")

# A parcellation: ~95 contiguous boxes clipped to the brain.
lab = np.zeros(brain.shape, dtype=np.int16)
idx = np.argwhere(brain)
key = (idx[:, 0] // 12) * 100 + (idx[:, 1] // 14) * 10 + (idx[:, 2] // 11)
u, inv = np.unique(key, return_inverse=True)
lab[brain] = inv + 1
nib.save(nib.Nifti1Image(lab, AFF, dtype=np.int16), f"{B}/labels.nii")

# The same data under an axis-permuted affine, for the reorientation probe.
AFF_ODD = np.array([[0., 0, 2, -50], [2., 0, 0, -60], [0, -2., 0, 70], [0, 0, 0, 1]])
nib.save(nib.Nifti1Image(d, AFF_ODD), f"{B}/odd.nii")

# Three sizes for the connected-component scaling probe.
rng2 = np.random.default_rng(3)
for name, sh in [("small", (40, 48, 38)), ("typical", (91, 109, 91)),
                 ("hires", (182, 218, 182))]:
    nib.save(nib.Nifti1Image(rng2.random(sh).astype(np.float32), np.diag([2., 2, 2, 1])),
             f"{B}/scale_{name}.nii")

# Sphere centres inside the brain, written in both conventions.
rng3 = np.random.default_rng(9)
pick = idx[rng3.choice(len(idx), 200, replace=False)]
np.savetxt(f"{B}/sphere_centres_world.txt", nib.affines.apply_affine(AFF, pick))
np.savetxt(f"{B}/sphere_centres_vox.txt", pick + 1, fmt="%d")   # 1-based for R

print(f"brain voxels {brain.sum()}  parcels {len(u)}  shape {SHAPE}  timepoints {T}")
