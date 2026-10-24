#!/usr/bin/env python3
"""Markdown yield table for all 2017 control-region categories.

Reads the per-variable yield logs that myControlPlots.C writes next to each plot and uses the
pT_MET one (pT_MET >= 0, so the histogram underflow is empty; overflow is included by the macro).
Usage: python yield_table.py [plots_dir] > yields_2017.md
"""
import glob
import os
import sys

PLOTS = sys.argv[1] if len(sys.argv) > 1 else "/depot/cms/users/shar1172/HZZ2l2nu_skim/controlplots_2017_AN325/plots"
GROUPS = ["Z+jets", "TTbar", "Diboson", "SingleTop", "TTV", "Triboson", "W+jets"]


def read(log):
    y = {}
    for line in open(log):
        if "\t=" in line:
            k, v = line.split("\t=")
            y[k.strip()] = float(v)
    return y


print("| category | data | total MC | data/MC | " + " | ".join(GROUPS) + " |")
print("|---|---:|---:|---:|" + "---:|" * len(GROUPS))
yields = {}
for region in ("LowMET", "mllSideband", "emu", "alphaSB"):
    for d in sorted(glob.glob(os.path.join(PLOTS, region, "*"))):
        tag = os.path.basename(d)
        logs = glob.glob(os.path.join(d, "*_pT_MET.log"))
        if not logs:
            continue
        y = read(logs[0])
        yields[tag] = y
        row = [tag, "%.0f" % y.get("data", 0), "%.1f" % y.get("Total MC", 0), "%.3f" % y.get("data/MC", 0)]
        row += ["%.1f" % y.get(g, 0) for g in GROUPS]
        print("| " + " | ".join(row) + " |")

# alpha method (AN-2016/325 Eq. 7-8): alpha_ll = N_ll / N_emu in the m_ll sidebands with MET > 70 and >= 1 b-tag;
# predicted non-resonant ee / mumu in the signal region = alpha_ll * N_emu (emu, on-Z, full selection, MET > 125).
sb = {f: yields.get("alphaSB_incl_" + f) for f in ("ee", "mumu", "emu")}
if all(sb.values()) and sb["emu"].get("data", 0) > 0:
    n_in = sum(yields[t].get("data", 0) for t in yields if t.startswith("emu_"))
    print("\n| alpha (AN-2016/325 5.1) | data | MC | AN-2016/325 data |")
    print("|---|---:|---:|---:|")
    for f, an in (("ee", "0.369 +- 0.006"), ("mumu", "0.683 +- 0.009")):
        nd, ne = sb[f].get("data", 0), sb["emu"].get("data", 0)
        md, me = sb[f].get("Total MC", 0), sb["emu"].get("Total MC", 0)
        err = (nd / ne) * (1.0 / nd + 1.0 / ne) ** 0.5 if nd > 0 else 0
        print("| alpha_%s | %.3f +- %.3f | %.3f | %s |" % (f, nd / ne, err, md / me if me else 0, an))
    print("\nN_emu^in (data, emu on-Z, MET > 125, all jet categories) = %.0f" % n_in)
