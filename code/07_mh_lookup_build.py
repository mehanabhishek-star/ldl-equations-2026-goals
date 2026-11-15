"""
07_mh_lookup_build.py
Exact Martin-Hopkins lookup for bridged triglycerides, recovered from NCHS's own values.

NHANES computes Martin-Hopkins LDL-C (LBDLDLM) as non-HDL-C - TG/factor, with the factor
taken from a 360-cell table (60 TG strata x 6 non-HDL-C strata) that the NCHS documentation
does not reproduce.
Because TG, non-HDL-C and LBDLDLM are integers, the VLDL-C estimate v = non-HDL-C - LBDLDLM
is the SAME integer for every participant with the same integer TG in the same non-HDL-C
band (<100, 100-129, 130-159, 160-189, 190-219, >=220). This script verifies that
(0 conflicts required) and writes v for every (band, TG) cell, so that Martin-Hopkins LDL-C
at a bridged TG is non-HDL-C - v(band, bridged TG): the value NCHS itself would report.

Cells with no NHANES participant are filled from the nearest observed TG values in the same
band (equal neighbours -> that value; otherwise interpolated factor) and flagged exact = 0.

Inputs : NHANES 2017-March 2020 (P_TRIGLY, P_TCHOL, P_HDL) and 2021-2023 (TRIGLY_L, TCHOL_L, HDL_L)
Output : results/mh_vldl_lookup.csv   (band 0-5, tg 10-399, v, exact, n_donors)
"""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
L = ROOT / "data" / "nhanes_2021_2023"
P = ROOT / "data" / "nhanes_2017_2020"
OUT = ROOT / "results" / "mh_vldl_lookup.csv"

def load(tr, tc, hd, tgvar):
    a = pd.read_sas(tr, format="xport")
    b = pd.read_sas(tc, format="xport")
    c = pd.read_sas(hd, format="xport")
    d = a.merge(b[["SEQN", "LBXTC"]], on="SEQN").merge(c[["SEQN", "LBDHDD"]], on="SEQN")
    return d.rename(columns={tgvar: "TG"})[["SEQN", "TG", "LBDLDLM", "LBXTC", "LBDHDD"]]

d = pd.concat([load(P / "P_TRIGLY.xpt", P / "P_TCHOL.xpt", P / "P_HDL.xpt", "LBXTR").assign(cyc=1),
               load(L / "TRIGLY_L.xpt", L / "TCHOL_L.xpt", L / "HDL_L.xpt", "LBXTLG").assign(cyc=2)])
d = d.dropna(subset=["TG", "LBDLDLM", "LBXTC", "LBDHDD"])
d["NH"] = d.LBXTC - d.LBDHDD
d["v"] = d.NH - d.LBDLDLM
assert ((d.TG % 1 == 0) & (d.NH % 1 == 0) & (d.LBDLDLM % 1 == 0)).all()
d["band"] = pd.cut(d.NH, [-1, 99, 129, 159, 189, 219, 9999], labels=False)

g = d.groupby(["band", "TG"]).v.agg(["nunique", "size", "first"])
assert (g["nunique"] == 1).all(), "v is not unique within (band, TG)"
both = d.groupby(["band", "TG", "cyc"]).v.first().unstack().dropna()
assert (both[1] == both[2]).all(), "cycles disagree"
print(f"records {len(d):,}; observed cells {len(g):,}; conflicts 0; cells in both cycles {len(both)} (all agree)")

obs = {(int(b), int(t)): (int(r["first"]), int(r["size"])) for (b, t), r in g.iterrows()}
rows = []
for b in range(6):
    ts = sorted(t for (bb, t) in obs if bb == b)
    for t in range(10, 400):
        if (b, t) in obs:
            rows.append((b, t, obs[(b, t)][0], 1, obs[(b, t)][1]))
            continue
        lo = [x for x in ts if x < t]
        hi = [x for x in ts if x > t]
        if lo and hi:
            t1, t2 = lo[-1], hi[0]
            v1, v2 = obs[(b, t1)][0], obs[(b, t2)][0]
            if v1 == v2:
                v = v1
            else:
                f = (t1 / v1) + (t2 / v2 - t1 / v1) * (t - t1) / (t2 - t1)
                v = int(np.floor(t / f + 0.5))
        elif lo:
            t1 = lo[-1]
            v = int(np.floor(t / (t1 / obs[(b, t1)][0]) + 0.5))
        else:
            t2 = hi[0]
            v = int(np.floor(t / (t2 / obs[(b, t2)][0]) + 0.5))
        rows.append((b, t, v, 0, 0))
out = pd.DataFrame(rows, columns=["band", "tg", "v", "exact", "n_donors"])
out.to_csv(OUT, index=False)
print("cells written", len(out), "| exact", int(out.exact.sum()), "| filled", int((out.exact == 0).sum()))
print(out.groupby("band").exact.mean().round(3).to_dict())
