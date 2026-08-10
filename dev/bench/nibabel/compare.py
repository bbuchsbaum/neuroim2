"""Compare neuroim2's probe results against nibabel's.

Run order:
    python3 make_probes.py                 # writes probes/ + probes/expected.json
    Rscript  probe_r.R > probes/r_results.json
    python3 compare.py

Only order-invariant summaries (min/max/sum/counts) and element [0] are
compared: nibabel flattens C-order and R flattens Fortran-order, so
"the second element" is a different voxel in the two languages and is not
a meaningful comparison.
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, "probes")

exp = json.load(open(os.path.join(P, "expected.json")))
got = {r["name"]: r for r in json.load(open(os.path.join(P, "r_results.json")))}


def close(a, b, tol=1e-6):
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a, float) and math.isnan(a):
        return isinstance(b, float) and math.isnan(b)
    return abs(a - b) <= tol * max(1.0, abs(b))


rows = []
for k in sorted(exp):
    v = exp[k]
    base = os.path.basename(v["path"])
    for suf in (".nii.gz", ".nii", ".hdr"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    g = got.get(base)
    if "error" in v:
        rows.append((k, "nibabel-cannot-load", v["error"][:80]))
        continue
    if g is None:
        rows.append((k, "NO-RESULT", ""))
        continue
    if not g["ok"]:
        rows.append((k, "neuroim2 ERROR", g["msg"][:95]))
        continue
    d = v["data"]
    if d["kind"] != "numeric":
        rows.append((k, "read w/o error", f"(reference is {d['kind']}; not compared)"))
        continue
    notes = []
    if list(g["shape"]) != list(d["shape"]):
        notes.append(f"shape {g['shape']} vs {d['shape']}")
    for f in ("first", "min", "max", "sum"):
        if not close(g[f], d[f]):
            notes.append(f"{f} {g[f]} vs {d[f]}")
    for f in ("n_nan", "n_inf"):
        if g[f] != d[f]:
            notes.append(f"{f} {g[f]} vs {d[f]}")
    aff_e = sum(v["affine"], [])
    if not all(abs(x - y) < 1e-4 for x, y in zip(g["affine"], aff_e)):
        notes.append("AFFINE differs")
    rows.append((k, "match" if not notes else "MISMATCH", "; ".join(notes)[:120]))

w = max(len(r[0]) for r in rows)
print(f"{'probe':<{w}}  {'result':<22} detail")
print("-" * (w + 24 + 40))
for k, s, n in rows:
    print(f"{k:<{w}}  {s:<22} {n}")
n_match = sum(1 for r in rows if r[1] == "match")
print(f"\n{n_match}/{len(rows)} probes agree with nibabel")
sys.exit(0 if n_match == len(rows) else 1)
