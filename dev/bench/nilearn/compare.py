"""Diff neuroim2's results against nilearn's and print the report.

    NEUROIM2_NL_DIR=/tmp/nl python3 compare.py
"""
import json, os
import numpy as np

B = os.environ.get("NEUROIM2_NL_DIR", "/tmp/nl")
ref = json.load(open(f"{B}/ref.json"))
got = json.load(open(f"{B}/res.json"))


def head(t):
    print(f"\n{t}\n" + "-" * len(t))


def ratio(a, b):
    return "n/a" if not b else f"{a / b:.0f}x" if a / b >= 10 else f"{a / b:.1f}x"


head("Smoothing: what came out, not what was asked for")
print(f"{'call':<44}{'PSF FWHM':>11}{'requested':>11}{'shortfall':>11}")
for f, v in ref["psf"].items():
    print(f"{'nilearn smooth_img(fwhm=' + f + ')':<44}{v:>9.2f}mm{float(f):>9.2f}mm{0:>10.0f}%")
for k, v in got["psf"].items():
    sg, w = k.replace("sigma", "").split("_window")
    lab = f"neuroim2 gaussian_blur(sigma={sg}, window={w})"
    print(f"{lab:<44}{v['fwhm']:>9.2f}mm{v['requested']:>9.2f}mm"
          f"{100 * (1 - v['fwhm'] / v['requested']):>10.0f}%")
print("\n  Timing. window=2 is the usual call but under-smooths, so it is doing less")
print("  work than nilearn; window=4 is the matched-kernel comparison.")
for lab, k in (("3-D", "t_smooth_3d"), ("4-D", "t_smooth_4d")):
    print(f"  {lab} smooth  neuroim2 window=2 {got[k]:.3f}s  window=4 {got[k + '_matched']:.3f}s"
          f"   nilearn {ref[k]:.3f}s   ({ratio(got[k + '_matched'], ref[k])} at matched kernel)")

head("Resampling")
print(f"  identity resample max|diff|: neuroim2 {got['identity_resample_maxdiff']:.3e}"
      f"   nilearn {ref['identity_resample_maxdiff']:.3e}")
for kind, v in got["resample"].items():
    print(f"  2mm->3mm {kind:<7} vs nilearn: cor {v['cor']:.6f}  RMS {v['rms']:.4f}"
          f"  max|diff| {v['max_abs']:.3f}")
print(f"  nearest on labels: {got['nearest_labels_in']} in -> neuroim2 {got['nearest_labels_out']},"
      f" nilearn {ref['nearest_labels_out']}")
print(f"  timing: neuroim2 {got['t_resample']:.3f}s  nilearn {ref['t_resample']:.3f}s"
      f"  ({ratio(got['t_resample'], ref['t_resample'])})")
nt = got["naive_target"]
print(f"\n  NeuroSpace(dim, spacing, origin) as a resample target:")
print(f"    source mean {nt['source_mean']:.4f} -> resampled mean {nt['resampled_mean']:.4f}"
      f"   {'<- target grid barely overlaps the source, silently' if nt['resampled_mean'] < 0.5 * nt['source_mean'] else ''}")

head("Reorientation to canonical (RAS)")
for who, v in (("nilearn reorder_img", ref["reorder"]), ("neuroim2 as_canonical", got["reorder"])):
    print(f"  {who:<24} {v['in_axcodes']} {v['in_shape']} -> {v['out_axcodes']} {v['out_shape']}"
          f"   nonzero {v['nonzero_in']} -> {v['nonzero_out']}"
          f" ({100 * v['nonzero_out'] / v['nonzero_in']:.1f}% kept)"
          f"   pure permutation: {v['values_are_a_permutation']}   {v['seconds']:.3f}s")

head("Connected components")
print(f"  {int(got['cc_in_mask'])} voxels above threshold")
for k in ref["cc"]:
    r, g = ref["cc"][k], got["cc"][k]
    agree = "match" if (r["n"] == g["n"] and r["largest"] == g["largest"]) else "MISMATCH"
    print(f"  {k:<12} scipy {r['n']:>6} comps / largest {r['largest']:>7} / {r['seconds']:7.3f}s"
          f"   neuroim2 {int(g['n']):>6} / {int(g['largest']):>7} / {g['seconds']:8.3f}s"
          f"   {agree}  ({ratio(g['seconds'], r['seconds'])})")
if "cc_scaling" in ref and "cc_scaling" in got:
    print(f"\n  {'size':<10}{'in-mask':>10}{'comps':>9}{'scipy':>10}{'neuroim2':>12}{'ratio':>9}")
    for k in ref["cc_scaling"]:
        r, g = ref["cc_scaling"][k], got["cc_scaling"][k]
        ok = "" if r["n"] == g["n"] else "  MISMATCH"
        print(f"  {k:<10}{r['in_mask']:>10}{r['n']:>9}{r['seconds']:>9.3f}s"
              f"{g['seconds']:>11.1f}s{ratio(g['seconds'], r['seconds']):>9}{ok}")

head("Brain masking from the data")
r, g = ref["epi_mask"], got["epi_mask"]
print(f"  nilearn compute_epi_mask  {r['voxels']:>7} voxels  Dice {r['dice']:.4f}  {r['seconds']:7.3f}s")
print(f"  neuroim2 automask         {int(g['voxels']):>7} voxels  Dice {g['dice']:.4f}  {g['seconds']:7.3f}s"
      f"   ({ratio(g['seconds'], r['seconds'])})")

head("Signal extraction (where neuroim2 is the faster side)")
for rad in got["spheres"]:
    r, g = ref["spheres"][rad], got["spheres"][rad]
    py = np.load(f"{B}/py_spheres_{rad}.npy")
    rr = np.loadtxt(f"{B}/r_spheres_{rad}.txt")
    print(f"  spheres r={rad}mm   nilearn {r['seconds']:6.3f}s   neuroim2 {g['seconds']:6.3f}s"
          f"   ({ratio(r['seconds'], g['seconds'])} faster)   max|diff| {np.abs(py - rr).max():.2e}")
py = np.load(f"{B}/py_labels.npy"); rr = np.loadtxt(f"{B}/r_labels.txt")
print(f"  parcel means      nilearn {ref['labels_signal']['seconds']:6.3f}s"
      f"   neuroim2 {got['labels_signal']['seconds']:6.3f}s   max|diff| {np.abs(py - rr).max():.2e}")
print(f"  masked extract    nilearn {ref['masked_extract']['seconds']:6.3f}s"
      f"   neuroim2 {got['masked_extract']['seconds']:6.3f}s")
py = np.load(f"{B}/py_mean.npy").ravel(order="F"); rr = np.loadtxt(f"{B}/r_mean.txt")
print(f"  mean image        max|diff| {np.abs(py - rr).max():.2e}")

head("Searchlight neighbourhoods")
for rad in got["searchlight"]:
    r, g = ref["searchlight"][rad], got["searchlight"][rad]
    print(f"  r={rad}mm  sklearn {r['centres']} centres mean {r['mean']:.1f} median {r['median']}"
          f" / {r['seconds']:.2f}s     neuroim2 {int(g['centres'])} mean {g['mean']:.1f}"
          f" median {int(g['median'])} / {g['seconds']:.2f}s")
cs = got["centre_set"]
print(f"\n  centre sets: searchlight() {int(cs['searchlight'])}, "
      f"searchlight_coords() {int(cs['searchlight_coords'])}, in-mask voxels {int(cs['in_mask'])}"
      + ("   <- the two iterators disagree"
         if cs["searchlight"] != cs["searchlight_coords"] else ""))

head("Utilities (seconds)")
print(f"  {'operation':<20}{'neuroim2':>10}{'nilearn':>10}{'ratio':>9}")
for k in ref["util"]:
    print(f"  {k:<20}{got['util'][k]:>10.3f}{ref['util'][k]:>10.3f}"
          f"{ratio(got['util'][k], ref['util'][k]):>9}")
print()
