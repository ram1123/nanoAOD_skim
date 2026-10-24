#!/usr/bin/env python3
"""
Build per-Run-2-year DAS sample lists for the X/H->ZZ->2l2nu skim from the
sample list of CMS AN-2016/325 (Tables 1, 2, 3, 52), extended with the
backgrounds added for the UL analysis.

AN-2016/325 is a 2016 (pre-UltraLegacy) note.  This script maps every AN
process to its UltraLegacy equivalent and resolves the exact DAS dataset
name(s) with `dasgoclient`, per year:

    2016preVFP  -> RunIISummer20UL16NanoAODAPV<nv> / HIPM_UL2016_*  (eras B-F)
    2016postVFP -> RunIISummer20UL16NanoAOD<nv>    / UL2016_*       (eras F-H, no HIPM)
    2017        -> RunIISummer20UL17NanoAOD<nv>    / UL2017_*       (eras B-F)
    2018        -> RunIISummer20UL18NanoAOD<nv>    / UL2018_*       (eras A-D)

Writes (into this directory):
    sample_list_AN2016_325_<v9|v15>_<year>.dat        - MC + DATA (one Slurm submission)
    sample_list_AN2016_325_<v9|v15>_<year>_data.dat   - DATA only

Format: one DAS dataset per line; '#' in column 0 = comment / section header
(compatible with scripts/condor/condor_setup_lxplus.py and
scripts/slurm/slurm_setup.py).  A process with no DAS match is written as
`# NOT IN DAS: <pattern>`; signal systematic variants are written commented.

Dataset choice per primary dataset: the standard NanoAOD processing (JMENano /
BTVNano only if nothing else exists), the highest -vN, and every _extN
extension (independent statistics; normalise with the summed genEventSumw).

Requires: a valid VOMS proxy (`voms-proxy-init -voms cms`) and `dasgoclient`
(on lxplus / cvmfs).  Run:

    python3 input_data_files/make_AN2016_325_lists.py --nano v15
    python3 input_data_files/make_AN2016_325_lists.py --nano v15 --years 2018
    python3 input_data_files/make_AN2016_325_lists.py            # NanoAODv9
"""
import argparse
import fnmatch
import re
import os
import shutil
import subprocess
import sys

OUTDIR = os.path.dirname(os.path.abspath(__file__))
DASGOCLIENT = shutil.which("dasgoclient") or "/cvmfs/cms.cern.ch/common/dasgoclient"

YEARS = {
    "2016preVFP": dict(
        mc="RunIISummer20UL16NanoAODAPV{NV}",
        data_run="Run2016", data_yy="2016",
        data_keep="HIPM",          # eras B(ver1,ver2),C,D,E,F  (HIPM / APV)
        data_eras="BCDEF",
        tier_note="UL16 pre-VFP (APV / HIPM)",
    ),
    "2016postVFP": dict(
        mc="RunIISummer20UL16NanoAOD{NV}",
        data_run="Run2016", data_yy="2016",
        data_keep="NOT_HIPM",      # eras F,G,H
        data_eras="FGH",
        tier_note="UL16 post-VFP",
    ),
    "2017": dict(
        mc="RunIISummer20UL17NanoAOD{NV}",
        data_run="Run2017", data_yy="2017",
        data_keep="ALL",
        # Run2017G (5 TeV) / H (low pileup) are outside the UL2017 golden JSON
        # (Cert_294927-306462, last run 306460): every event would be rejected.
        data_eras="BCDEF",
        tier_note="UL17",
    ),
    "2018": dict(
        mc="RunIISummer20UL18NanoAOD{NV}",
        data_run="Run2018", data_yy="2018",
        data_keep="ALL",
        data_eras="ABCD",
        tier_note="UL18",
    ),
}

# conditions / re-processing substrings we never want (primary-dataset level)
BAD_SUBSTR = (
    "PUForMUOVal", "PUForTRK", "PUForNanoMuon", "PUFor", "Pilot",
    "FlatPU", "PU35", "PU25", "LowPU", "LensingPU", "_BS20", "PrivateMC",
    "EpsilonPU", "forPOG", "NanoMuon",
)
# special NanoAOD flavours (processed-dataset level): only if nothing else exists
SPECIAL_NANO = ("JMENano", "BTVNano")
# tune / model systematic variants - dropped when a nominal exists
SYST_SUBSTR = (
    "TuneCP5CR1", "TuneCP5CR2", "TuneCP5up", "TuneCP5down", "TuneCP5Up",
    "TuneCP5Down", "TuneCP5_erdON", "CP5TuneUp", "CP5TuneDown", "_erdON",
    "_mtop1", "hdamp", "GEN-", "minloHJJ", "BGenFilter", "TuneCH3", "herwig",
    "evtgen", "withDipoleRecoil",
)

# --------------------------------------------------------------------------
# AN-2016/325 sample list -> UL primary-dataset name patterns (fnmatch on /PD).
#   entry:  (label, [patterns], mode)
#   mode "first"  : alternative names for ONE sample, first hit wins
#   mode "all"    : each pattern is a distinct sample (HT/pT bins)
#   mode "signal" : every mass point found; systematic variants written commented
# tune is TuneCP5 in UL (the AN used TuneCUETP8M1 / RunIISummer16).
# --------------------------------------------------------------------------
def _bins(stem, bins, tail):
    return ["/%s%s%s" % (stem, b, tail) for b in bins]


MC = [
    ("### W->lnu  (AN Table 2, 52)", None, None),
    ("WJetsToLNu inclusive (amcatnloFXFX)",
     ["/WJetsToLNu_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("WJetsToLNu HT-binned (madgraphMLM)",
     _bins("WJetsToLNu_HT-", ("70To100", "100To200", "200To400", "400To600", "600To800",
                              "800To1200", "1200To2500", "2500ToInf"),
           "_TuneCP5_13TeV-madgraphMLM-pythia8"), "all"),

    ("### Z->ll  Drell-Yan  (AN Table 2)", None, None),
    ("DYJetsToLL M-10to50 (amcatnloFXFX)",
     ["/DYJetsToLL_M-10to50_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("DYJetsToLL M-50 (amcatnloFXFX)",
     ["/DYJetsToLL_M-50_TuneCP5*_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("DYJetsToLL LHEFilterPtZ-binned (amcatnloFXFX)",
     _bins("DYJetsToLL_LHEFilterPtZ-", ("0To50", "50To100", "100To250", "250To400",
                                        "400To650", "650ToInf"),
           "_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8"), "all"),

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
     ["/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"], "first"),
    ("ST tW antitop 5f inclusiveDecays (powheg)",
     ["/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"], "first"),
    ("ST t-channel top 4f InclusiveDecays (powheg-madspin)",
     ["/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"], "first"),
    ("ST t-channel antitop 4f InclusiveDecays (powheg-madspin)",
     ["/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"], "first"),
    ("ST s-channel 4f leptonDecays (amcatnlo)",
     ["/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8"], "first"),
    ("ST s-channel 4f hadronicDecays (amcatnlo)",
     ["/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8"], "first"),
    ("ST tW top/antitop 5f NoFullyHadronicDecays (powheg)",
     ["/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
      "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"], "all"),
    ("tZq ll 4f ckm NLO (amcatnlo)",
     ["/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8"], "first"),

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
    ("VBFToContinToZZTo2l2nu (JHUGen+MCFM701)  -- EW qqZZjj continuum, VBF signal interference",
     ["/VBFToContinToZZTo2l2nu_TuneCP5*_13TeV-jhugenv730-mcfm701-pythia8"], "first"),
    ("GluGluToWW -> 2l2nu (MCFM701)  -- gg->WW continuum",
     _bins("GluGluToWWTo", ("ENEN", "MNMN", "ENMN"), "_TuneCP5_13TeV_MCFM701_pythia8"), "all"),
    ("WZTo1L1Nu2Q (amcatnloFXFX)",
     ["/WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("WZTo1L3Nu (amcatnloFXFX)",
     ["/WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("WZTo2Q2Nu (amcatnloFXFX-madspin)",
     ["/WZTo2Q2Nu_4f_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"], "first"),
    ("ZZTo2Nu2Q (amcatnloFXFX)",
     ["/ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8"], "first"),
    ("ZZTo4L (powheg)",
     ["/ZZTo4L_TuneCP5_13TeV_powheg_pythia8"], "first"),

    ("### Tribosons  (AN Table 2, 52)", None, None),
    ("ZZZ (amcatnlo)", ["/ZZZ_TuneCP5*_13TeV-amcatnlo-pythia8"], "first"),
    ("WZZ (amcatnlo)", ["/WZZ_TuneCP5*_13TeV-amcatnlo-pythia8"], "first"),
    ("WWZ (amcatnlo)",
     ["/WWZ_4F_TuneCP5*_13TeV-amcatnlo-pythia8",
      "/WWZ_TuneCP5*_13TeV-amcatnlo-pythia8"], "first"),
    ("WWW (amcatnlo)",
     ["/WWW_4F_TuneCP5*_13TeV-amcatnlo-pythia8"], "first"),

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

    ("### Signal  (AN section 3 - POWHEG + JHUGen, high-mass grid; every mass point in DAS)",
     None, None),
    ("GluGluHToZZTo2L2Nu M*",
     ["/GluGluHToZZTo2L2Nu_M*_13TeV_powheg2*_pythia8"], "signal"),
    ("VBF HToZZTo2L2Nu M*",
     ["/VBF_HToZZTo2L2Nu_M*_13TeV_powheg2*_pythia8"], "signal"),
]

DATA_PD_RUN2 = [("DoubleMuon", ""), ("DoubleEG", ""), ("SingleMuon", ""),
                ("SingleElectron", ""), ("MuonEG", "  (emu control region)"),
                ("SinglePhoton", "  (HLT_Photon* in Triggers_HZZ2l2nu_SR for high-pT electrons; gamma CR)")]
DATA_PD_2018 = [("DoubleMuon", ""), ("EGamma", "  (DoubleEG + SingleElectron + SinglePhoton)"),
                ("SingleMuon", ""), ("MuonEG", "  (emu control region)")]


def das(query):
    try:
        out = subprocess.run([DASGOCLIENT, "--query", query],
                             capture_output=True, text=True, timeout=600)
    except Exception as e:                                    # noqa: BLE001
        sys.stderr.write("dasgoclient error for %r: %s\n" % (query, e))
        return []
    if out.returncode != 0:
        sys.stderr.write("dasgoclient failed for %r: %s\n" % (query, out.stderr.strip()))
    return [l.strip() for l in out.stdout.splitlines() if l.strip().startswith("/")]


def _pd(d):
    return d.split("/")[1]


def _proc(d):
    return d.split("/")[2]


def select_versions(dsets):
    """For one primary dataset: standard NanoAOD flavour, then per extension
    (none, _ext1, ...) the plain processing string and the highest -vN."""
    std = [d for d in dsets if not any(s in _proc(d) for s in SPECIAL_NANO)] or dsets
    groups = {}
    for d in std:
        m = re.match(r'(.*)-v(\d+)$', _proc(d))
        stem, ver = (m.group(1), int(m.group(2))) if m else (_proc(d), 0)
        ext = re.search(r'_ext\d+', stem)
        ext = ext.group(0) if ext else ""
        # plain processing ("<campaign>-150X_...") beats prefixed ones ("4cores5k_150X_...")
        plain = 1 if re.search(r'NanoAOD(APV)?v\d+-\d+X_', stem) else 0
        key = (plain, ver)
        if ext not in groups or key > groups[ext][0]:
            groups[ext] = (key, d)
    return [groups[e][1] for e in sorted(groups)]


def match_pds(catalog, pattern):
    pds = sorted({_pd(d) for d in catalog if fnmatch.fnmatchcase("/" + _pd(d), pattern)})
    pds = [p for p in pds if not any(b in p for b in BAD_SUBSTR)] or pds
    return pds


def pick(catalog, pds):
    """drop tune/model systematic variants and PSWeights duplicates when a
    nominal exists; return the selected datasets of the remaining PDs."""
    if not pds:
        return []
    nominal = [p for p in pds if not any(s in p for s in SYST_SUBSTR)] or pds
    no_ps = [p for p in nominal if "PSWeights" not in p] or nominal
    out = []
    for p in no_ps:
        out += select_versions([d for d in catalog if _pd(d) == p])
    return out


def _mass(pd):
    m = re.search(r'_M-?(\d+)', pd)
    return int(m.group(1)) if m else 0


def build_mc(year, nv):
    cfg = YEARS[year]
    campaign = cfg["mc"].format(NV=nv)
    catalog = das("dataset=/*/%s-*/NANOAODSIM" % campaign)
    if not catalog:
        sys.exit("No %s datasets returned by DAS (proxy / network?)" % campaign)
    L = []
    ok = bad = 0
    for label, patterns, mode in MC:
        if patterns is None:
            L += ["#", "#### " + label[4:]]
            continue
        L.append("# %s" % label)
        if mode == "signal":
            pds = match_pds(catalog, patterns[0])
            nom = sorted([p for p in pds if not any(s in p for s in SYST_SUBSTR)], key=_mass)
            syst = sorted([p for p in pds if any(s in p for s in SYST_SUBSTR)], key=lambda p: (_mass(p), p))
            for p in nom:
                L += select_versions([d for d in catalog if _pd(d) == p])
            L.append("#   masses in DAS: %s" % " ".join(str(_mass(p)) for p in nom))
            if syst:
                L.append("#   systematic / alternative samples (tune Up/Down, minloHJJ, dipole recoil): uncomment to include")
                for p in syst:
                    L += ["#" + d for d in select_versions([d for d in catalog if _pd(d) == p])]
            if nom:
                ok += 1
            else:
                L.append("# NOT IN DAS: dataset=%s/%s-*/NANOAODSIM" % (patterns[0], campaign))
                bad += 1
        elif mode == "all":
            for pat in patterns:
                hits = pick(catalog, match_pds(catalog, pat))
                if hits:
                    L += hits
                    ok += 1
                else:
                    L.append("# NOT IN DAS: dataset=%s/%s-*/NANOAODSIM" % (pat, campaign))
                    bad += 1
        else:
            hits = None
            for pat in patterns:
                hits = pick(catalog, match_pds(catalog, pat))
                if hits:
                    break
            if hits:
                L += hits
                ok += 1
            else:
                L.append("# NOT IN DAS: dataset=%s/%s-*/NANOAODSIM" % (patterns[-1], campaign))
                bad += 1
    return L, ok, bad


def build_data(year, nv):
    cfg = YEARS[year]
    yy = cfg["data_yy"]
    mid = "MiniAODv2_NanoAODv9" if nv == "v9" else "NanoAODv15"
    keep = cfg["data_keep"]
    pds = DATA_PD_2018 if year == "2018" else DATA_PD_RUN2
    catalog = das("dataset=/*/%s*UL%s_%s*/NANOAOD" % (cfg["data_run"], yy, mid))
    L = []
    if year == "2017":
        L += ["# Run2017B-F only: Run2017G (5 TeV) / H (low pileup) are outside the UL2017 golden JSON",
              "# (Cert_294927-306462, last run 306460), so every event would be rejected."]
    ok = bad = 0
    for pd, note in pds:
        best = {}
        for d in catalog:
            if _pd(d) != pd:
                continue
            m = re.match(r'%s(\w)-(.*)-v(\d+)$' % cfg["data_run"], _proc(d))
            if not m or m.group(1) not in cfg["data_eras"]:
                continue
            hipm = "HIPM" in m.group(2)
            if (keep == "HIPM" and not hipm) or (keep == "NOT_HIPM" and hipm):
                continue
            # one dataset per era and processing string (2016B ver1 / ver2 are both kept)
            key, ver = (m.group(1), m.group(2)), int(m.group(3))
            if key not in best or ver > best[key][0]:
                best[key] = (ver, d)
        hits = [best[k][1] for k in sorted(best)]
        L.append("# %s%s" % (pd, note))
        if hits:
            L += hits
            ok += 1
        else:
            L.append("# NOT IN DAS: dataset=/%s/%s*UL%s_%s*/NANOAOD (eras %s%s)" % (
                pd, cfg["data_run"], yy, mid, cfg["data_eras"],
                ", %s filter" % keep if keep != "ALL" else ""))
            bad += 1
    return L, ok, bad


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
        cfg = YEARS[year]
        mc, a1, b1 = build_mc(year, args.nano)
        data, a2, b2 = build_data(year, args.nano)
        head = ["# X/H->ZZ->2l2nu  %%s sample list  -  %s  (%s, NanoAOD %s)" % (
                    year, cfg["tier_note"], args.nano),
                "# from CMS AN-2016/325 Tables 1, 2, 3, 52  ->  UltraLegacy equivalents (+ UL additions)",
                "# generated by input_data_files/make_AN2016_325_lists.py (re-run to refresh; edits are overwritten)",
                "#"]
        p1 = os.path.join(OUTDIR, "sample_list_AN2016_325_%s_%s.dat" % (args.nano, year))
        p2 = os.path.join(OUTDIR, "sample_list_AN2016_325_%s_%s_data.dat" % (args.nano, year))
        body = ([head[0] % "MC + DATA"] + head[1:] +
                ["# One list for one Slurm submission: post_proc.py sets isMC / golden JSON per file (\"/data/\" in LFN).",
                 "#"] + mc +
                ["#", "#### DATA  (AN-2016/325 Table 1 primary datasets -> UL%s NanoAOD %s)" % (
                    cfg["data_yy"], args.nano)] + data)
        open(p1, "w").write("\n".join(body) + "\n")
        open(p2, "w").write("\n".join([head[0] % "DATA"] + head[1:] + data) + "\n")
        print("%-12s  MC   %-40s  %3d ok / %2d not in DAS" % (year, os.path.basename(p1), a1, b1))
        print("%-12s  DATA %-40s  %3d ok / %2d not in DAS" % ("", os.path.basename(p2), a2, b2))


if __name__ == "__main__":
    main()
