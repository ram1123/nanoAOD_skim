#!/usr/bin/env python3
"""
Build per-Run-2-year DAS sample lists for the X/H->ZZ->2l2nu skim from the
sample list of CMS AN-2016/325 (Tables 1, 2, 3, 52).

AN-2016/325 is a 2016 (pre-UltraLegacy) note.  This script maps every AN
process to its **UltraLegacy NanoAODv9** equivalent and resolves the exact
DAS dataset name(s) with `dasgoclient`, per year:

    2016preVFP  -> RunIISummer20UL16NanoAODAPVv9 / HIPM_UL2016_MiniAODv2_NanoAODv9
    2016postVFP -> RunIISummer20UL16NanoAODv9    / UL2016_MiniAODv2_NanoAODv9 (no HIPM)
    2017        -> RunIISummer20UL17NanoAODv9    / UL2017_MiniAODv2_NanoAODv9
    2018        -> RunIISummer20UL18NanoAODv9    / UL2018_MiniAODv2_NanoAODv9

Writes (into this directory):
    sample_list_AN2016_325_<v9|v15>_<year>.dat        - MC, grouped by AN category
    sample_list_AN2016_325_<v9|v15>_<year>_data.dat   - Data primary datasets

Format: one DAS dataset per line; '#' in column 0 = comment / section header
(compatible with scripts/condor/condor_setup_lxplus.py and
scripts/slurm/slurm_setup.py).  Patterns with no DAS match are written as
`# UNRESOLVED: <pattern>` for the user to fix by hand.

Requires: a valid VOMS proxy (`voms-proxy-init -voms cms`) and `dasgoclient`
(on lxplus / cvmfs).  Run:

    python3 input_data_files/make_AN2016_325_lists.py            # NanoAODv9 (default)
    python3 input_data_files/make_AN2016_325_lists.py --nano v15 # if v15 re-nano wanted
    python3 input_data_files/make_AN2016_325_lists.py --years 2018
"""
import argparse
import re
import os
import shutil
import subprocess
import sys

OUTDIR = os.path.dirname(os.path.abspath(__file__))
DASGOCLIENT = shutil.which("dasgoclient") or "/cvmfs/cms.cern.ch/common/dasgoclient"

YEARS = {
    "2016preVFP": dict(
        mc="RunIISummer20UL16NanoAODAPV{NV}-*",
        data_run="Run2016", data_yy="2016",
        data_keep="HIPM",          # eras B(ver1,ver2),C,D,E,F  (HIPM / APV)
        tier_note="UL16 pre-VFP (APV / HIPM)",
    ),
    "2016postVFP": dict(
        mc="RunIISummer20UL16NanoAOD{NV}-*",
        data_run="Run2016", data_yy="2016",
        data_keep="NOT_HIPM",      # eras F,G,H
        tier_note="UL16 post-VFP",
    ),
    "2017": dict(
        mc="RunIISummer20UL17NanoAOD{NV}-*",
        data_run="Run2017", data_yy="2017",
        data_keep="ALL",           # eras B..F
        tier_note="UL17",
    ),
    "2018": dict(
        mc="RunIISummer20UL18NanoAOD{NV}-*",
        data_run="Run2018", data_yy="2018",
        data_keep="ALL",           # eras A..D
        tier_note="UL18",
    ),
}

# conditions / re-processing substrings we never want
BAD_SUBSTR = (
    "JMENano", "BTVNano", "PUForMUOVal", "PUForTRK", "PUForNanoMuon", "PUFor", "Pilot",
    "FlatPU", "PU35", "PU25", "LowPU", "LensingPU", "_BS20", "PrivateMC",
    "EpsilonPU", "forPOG", "NanoMuon",
)
# tune / model systematic variants - dropped when a nominal exists
SYST_SUBSTR = (
    "TuneCP5CR1", "TuneCP5CR2", "TuneCP5up", "TuneCP5down", "TuneCP5Up",
    "TuneCP5Down", "TuneCP5_erdON", "CP5TuneUp", "CP5TuneDown", "_erdON",
    "_mtop1", "hdamp", "GEN-", "minloHJJ", "BGenFilter", "TuneCH3", "herwig",
    "evtgen",
)

# --------------------------------------------------------------------------
# AN-2016/325 sample list -> UL primary-dataset name patterns.
#   entry:  (label, [patterns], mode)
#   mode "first" (default) : alternative names for ONE sample, first hit wins
#   mode "all"             : each pattern is a distinct sample (HT/pT bins)
# tune is TuneCP5 in UL (the AN used TuneCUETP8M1 / RunIISummer16).
# --------------------------------------------------------------------------
def _bins(stem, bins, tail):
    return ["/%s%s%s" % (stem, b, tail) for b in bins]


MC = [
    ("### W->lnu  (AN Table 2, 52)", None, None),
    ("WJetsToLNu inclusive (amcatnloFXFX)",
     ["/WJetsToLNu_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("WJetsToLNu HT-binned (madgraphMLM)",
     _bins("WJetsToLNu_HT-", ("100To200", "200To400", "400To600", "600To800",
                              "800To1200", "1200To2500", "2500ToInf"),
           "_TuneCP5_13TeV-madgraphMLM-pythia8"), "all"),

    ("### Z->ll  Drell-Yan  (AN Table 2)", None, None),
    ("DYJetsToLL M-10to50 (amcatnloFXFX)",
     ["/DYJetsToLL_M-10to50_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("DYJetsToLL M-50 (amcatnloFXFX)",
     ["/DYJetsToLL_M-50_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),

    ("### Z->nunu  HT>100  (AN Table 52 - instrumental-MET closure)", None, None),
    ("ZJetsToNuNu HT-binned (madgraphMLM)",
     _bins("ZJetsToNuNu_HT-", ("100To200", "200To400", "400To600", "600To800",
                               "800To1200", "1200To2500", "2500ToInf"),
           "_TuneCP5_13TeV-madgraphMLM-pythia8"), "all"),

    ("### tt + X  (AN Table 2)", None, None),
    ("TTTo2L2Nu (powheg)", ["/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"], "first"),
    ("TTWJetsToLNu (amcatnloFXFX-madspin)",
     ["/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"], "first"),
    ("TTZToLLNuNu M-10 (amcatnlo)",
     ["/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8"], "first"),

    ("### Single top  (AN Table 2, 52)", None, None),
    ("ST tW top 5f inclusiveDecays (powheg)",
     ["/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
      "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"], "first"),
    ("ST tW antitop 5f inclusiveDecays (powheg)",
     ["/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
      "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"], "first"),
    ("ST t-channel top 4f InclusiveDecays (powheg-madspin)",
     ["/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"], "first"),
    ("ST t-channel antitop 4f InclusiveDecays (powheg-madspin)",
     ["/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"], "first"),
    ("ST s-channel 4f leptonDecays (amcatnlo)",
     ["/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8"], "first"),

    ("### Dibosons  (AN Table 2, 52)", None, None),
    ("WZTo3LNu (powheg / amcatnloFXFX)",
     ["/WZTo3LNu_TuneCP5_13TeV-powheg-pythia8",
      "/WZTo3LNu_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("WZTo2L2Q  (UL: WZTo2Q2L_mllmin4p0, amcatnloFXFX)",
     ["/WZTo2Q2L_mllmin4p0_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("WWTo2L2Nu (powheg)",
     ["/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"], "first"),
    ("WWToLNuQQ  (UL: WWTo1L1Nu2Q_4f, amcatnloFXFX)",
     ["/WWToLNuQQ_TuneCP5_13TeV-powheg-pythia8",
      "/WWTo1L1Nu2Q_4f_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("ZZTo2L2Nu (powheg)  -- key irreducible",
     ["/ZZTo2L2Nu_TuneCP5*_13TeV_powheg_pythia8",
      "/ZZTo2L2Nu_TuneCP5*_13TeV-powheg-pythia8"], "first"),
    ("ZZTo2L2Q  (UL: ZZTo2Q2L_mllmin4p0, amcatnloFXFX)",
     ["/ZZTo2Q2L_mllmin4p0_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("GluGluToContinToZZTo2mu2nu (MCFM701)  -- gg->ZZ continuum",
     ["/GluGluToContinToZZTo2mu2nu_TuneCP5*_13TeV-mcfm701-pythia8"], "first"),
    ("GluGluToContinToZZTo2e2nu (MCFM701)  -- gg->ZZ continuum",
     ["/GluGluToContinToZZTo2e2nu_TuneCP5*_13TeV-mcfm701-pythia8"], "first"),

    ("### Tribosons  (AN Table 2, 52)", None, None),
    ("ZZZ (amcatnlo)", ["/ZZZ_TuneCP5*_13TeV-amcatnlo-pythia8"], "first"),
    ("WZZ (amcatnlo)", ["/WZZ_TuneCP5*_13TeV-amcatnlo-pythia8"], "first"),
    ("WWZ (amcatnlo)",
     ["/WWZ_4F_TuneCP5*_13TeV-amcatnlo-pythia8",
      "/WWZ_TuneCP5*_13TeV-amcatnlo-pythia8"], "first"),

    ("### gamma+jets  (AN Table 3, 52 - photon control region)", None, None),
    ("GJets HT-binned (madgraphMLM)",
     _bins("GJets_HT-", ("40To100", "100To200", "200To400", "400To600", "600ToInf"),
           "_TuneCP5_13TeV-madgraphMLM-pythia8"), "all"),

    ("### Top + gamma  (AN Table 3, 52)", None, None),
    ("TTGJets (amcatnloFXFX-madspin)",
     ["/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"], "first"),
    ("TGJets (amcatnlo-madspin)",
     ["/TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8",
      "/TGJets_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8"], "first"),

    ("### Z gamma  (AN Table 3, 52)", None, None),
    ("ZNuNuGJets MonoPhoton PtG-40to130",
     ["/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCP5*_13TeV-*pythia8"], "first"),
    ("ZNuNuGJets MonoPhoton PtG-130",
     ["/ZNuNuGJets_MonoPhoton_PtG-130_TuneCP5*_13TeV-*pythia8"], "first"),
    ("ZGamma -> ll gamma  (UL: ZGToLLG_01J_5f, amcatnloFXFX)",
     ["/ZGToLLG_01J_5f_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),

    ("### W gamma  (AN Table 3, 52)", None, None),
    ("WGToLNuG (madgraphMLM / amcatnloFXFX)",
     ["/WGToLNuG_TuneCP5*_13TeV-madgraphMLM-pythia8",
      "/WGToLNuG_01J_5f_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),

    ("### QCD  HT>100  (AN Table 52)", None, None),
    ("QCD HT-binned (madgraphMLM)",
     _bins("QCD_HT", ("100to200", "200to300", "300to500", "500to700",
                      "700to1000", "1000to1500", "1500to2000", "2000toInf"),
           "_TuneCP5*_13TeV-madgraph*-pythia8"), "all"),

    ("### QCD EMEnriched  (AN Table 52)", None, None),
    ("QCD Pt EMEnriched (pythia8)",
     _bins("QCD_Pt-", ("15to20", "20to30", "30to50", "50to80", "80to120",
                       "120to170", "170to300", "300toInf"),
           "_EMEnriched_TuneCP5_13TeV-pythia8"), "all"),

    ("### QCD MuEnriched  (AN Table 52)  -- UL uses capital 'To'; top bin is Pt-1000", None, None),
    ("QCD Pt MuEnrichedPt5 (pythia8)",
     _bins("QCD_Pt-", ("15To20", "20To30", "30To50", "50To80", "80To120",
                       "120To170", "170To300", "300To470", "470To600",
                       "600To800", "800To1000", "1000"),
           "_MuEnrichedPt5_TuneCP5_13TeV-pythia8"), "all"),

    ("### Signal  (AN section 3 - POWHEG + JHUGen, high-mass grid; trim to the wanted points)",
     None, None),
    ("GluGluHToZZTo2L2Nu M*",
     ["/GluGluHToZZTo2L2Nu_M*_TuneCP5_13TeV_powheg2_JHUGen*_pythia8"], "first"),
    ("VBF HToZZTo2L2Nu M*",
     ["/VBF_HToZZTo2L2Nu_M*_TuneCP5_13TeV_powheg2_JHUGen*_pythia8"], "first"),
]

DATA_PD_RUN2 = ["DoubleMuon", "DoubleEG", "SingleMuon", "SingleElectron",
                "MuonEG", "SinglePhoton"]
DATA_PD_2018 = ["DoubleMuon", "EGamma", "SingleMuon", "MuonEG"]


def das(query):
    try:
        out = subprocess.run([DASGOCLIENT, "--query", query],
                             capture_output=True, text=True, timeout=120)
    except Exception as e:                                    # noqa: BLE001
        sys.stderr.write("dasgoclient error for %r: %s\n" % (query, e))
        return []
    return [l.strip() for l in out.stdout.splitlines() if l.strip().startswith("/")]


def pick(results):
    """standard re-processing(s) only: drop specialised campaigns, then drop
    tune/model systematic variants when a nominal one exists."""
    if not results:
        return []
    good = [r for r in results if not any(b in r for b in BAD_SUBSTR)] or results
    nominal = [r for r in good if not any(s in r for s in SYST_SUBSTR)]
    good = nominal or good
    no_ps = [r for r in good if "PSWeights" not in r]
    good = no_ps or good
    base = [r for r in good if "_ext" not in r]
    return sorted(set(base or good))


def resolve_first(patterns, campaign, tier):
    for pat in patterns:
        hits = pick(das("dataset=%s/%s/%s" % (pat, campaign, tier)))
        if hits:
            return hits
    return None


def build_mc(year, nv):
    cfg = YEARS[year]
    campaign = cfg["mc"].format(NV=nv)
    L = ["# X/H->ZZ->2l2nu  MC sample list  -  %s  (%s, NanoAOD %s)" % (
        year, cfg["tier_note"], nv),
        "# from CMS AN-2016/325 Tables 2, 3, 52  ->  UltraLegacy equivalents",
        "# generated by input_data_files/make_AN2016_325_lists.py (re-run to refresh)",
        "#"]
    ok = bad = 0
    for label, patterns, mode in MC:
        if patterns is None:
            L += ["#", "#### " + label[4:]]
            continue
        if mode == "all":
            L.append("# %s" % label)
            for pat in patterns:
                hits = resolve_first([pat], campaign, "NANOAODSIM")
                if hits:
                    L += hits
                    ok += 1
                else:
                    L.append("# UNRESOLVED: dataset=%s/%s/NANOAODSIM" % (pat, campaign))
                    bad += 1
        else:
            hits = resolve_first(patterns, campaign, "NANOAODSIM")
            L.append("# %s" % label)
            if hits:
                L += hits
                ok += 1
            else:
                L.append("# UNRESOLVED: dataset=%s/%s/NANOAODSIM" % (patterns[-1], campaign))
                bad += 1
    path = os.path.join(OUTDIR, "sample_list_AN2016_325_%s_%s.dat" % (nv, year))
    open(path, "w").write("\n".join(L) + "\n")
    return path, ok, bad


def _latest_version(dsets):
    # per era (dataset name up to the trailing -vN), keep the highest N
    best = {}
    for d in dsets:
        m = re.match(r'(.*)-v(\d+)/NANOAOD$', d)
        key, ver = (m.group(1), int(m.group(2))) if m else (d, -1)
        if key not in best or ver > best[key][0]:
            best[key] = (ver, d)
    return [v[1] for v in sorted(best.values(), key=lambda t: t[1])]


def build_data(year, nv):
    cfg = YEARS[year]
    yy = cfg["data_yy"]
    mid = "MiniAODv2_NanoAODv9" if nv == "v9" else "NanoAODv15"
    keep = cfg["data_keep"]
    pds = DATA_PD_2018 if year == "2018" else DATA_PD_RUN2
    L = ["# X/H->ZZ->2l2nu  DATA sample list  -  %s  (%s, NanoAOD %s)" % (
        year, cfg["tier_note"], nv),
        "# AN-2016/325 Table 1 primary datasets  ->  UltraLegacy",
        "# generated by input_data_files/make_AN2016_325_lists.py",
        "#"]
    ok = bad = 0
    for pd in pds:
        q = "dataset=/%s/%s*UL%s_%s*/NANOAOD" % (pd, cfg["data_run"], yy, mid)
        hits = [h for h in sorted(set(das(q)))
                if "BTVNano" not in h and "JMENano" not in h]
        hits = _latest_version(hits)
        if keep == "HIPM":
            hits = [h for h in hits if "HIPM" in h]
        elif keep == "NOT_HIPM":
            hits = [h for h in hits if "HIPM" not in h]
        L.append("# %s" % pd)
        if hits:
            L += hits
            ok += 1
        else:
            L.append("# UNRESOLVED: %s%s" % (
                q, "  (%s filter)" % keep if keep != "ALL" else ""))
            bad += 1
    path = os.path.join(OUTDIR, "sample_list_AN2016_325_%s_%s_data.dat" % (nv, year))
    open(path, "w").write("\n".join(L) + "\n")
    return path, ok, bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nano", choices=["v9", "v15"], default="v9",
                    help="NanoAOD re-processing (default v9 = full UL coverage)")
    ap.add_argument("--years", nargs="+", default=list(YEARS), choices=list(YEARS))
    args = ap.parse_args()

    if subprocess.run(["voms-proxy-info", "-exists"],
                      capture_output=True).returncode != 0:
        sys.exit("No valid VOMS proxy - run: voms-proxy-init -voms cms")

    for year in args.years:
        p1, a1, b1 = build_mc(year, args.nano)
        p2, a2, b2 = build_data(year, args.nano)
        print("%-12s  MC   %-38s  %3d ok / %2d unresolved" % (year, os.path.basename(p1), a1, b1))
        print("%-12s  DATA %-38s  %3d ok / %2d unresolved" % ("", os.path.basename(p2), a2, b2))


if __name__ == "__main__":
    main()
