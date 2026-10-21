#!/usr/bin/env python3
"""Build slim plotting inputs (control plots / shapes) from the Slurm skims.

The skim area is <skimdir>/<sample>/*_Skim.root, as written by scripts/slurm/slurm_setup.py.
For every sample directory it writes ONE slim ROOT file
(<outdir>/<sample>.root, tree "Events") holding only
  * events passing passZZ2l2nuSelection || passZZ2l2nu_emuCR_Selection
    (every dilepton control region requires one of the two), and
  * the branches listed in BRANCHES (cuts, weights, plotted variables).

MC:   genEventSumw is summed over the Runs trees of ALL skim files of the sample
      (not just the selected events) and written to <outdir>/sumw.json.
Data: all primary datasets are chained in priority order and duplicate
      (run, luminosityBlock, event) are dropped (first occurrence kept), same
      logic as removeDuplicates.C. One file: <outdir>/<data_file>.

--mode gamma: single-photon CR skim (--channels 1g) instead; loose filter
passGammaCR_Selection && gamma_trigWindowFired, data = SinglePhoton only.

Runs with ROOT from the AF global pixi env (no CMSSW needed):
  /work/pixi/global/.pixi/envs/default/bin/python scripts/analysis/prepare_inputs.py \\
      --skimdir DIR --outdir DIR [--mode dilepton|gamma] [--data_file data_2017.root] \\
      [--only NAME ...] [--threads N] [--force]
Re-running skips samples whose output already exists (use --force to redo).
"""
import argparse
import glob
import json
import os
import sys

import ROOT

# Data primary datasets, in duplicate-removal priority order.
DATA_PDS = ["DoubleMuon", "DoubleEG", "MuonEG", "SingleMuon", "SingleElectron", "SinglePhoton"]

LOOSE = "passZZ2l2nuSelection || passZZ2l2nu_emuCR_Selection"

# --mode gamma (single-photon CR skim, --channels 1g)
GAMMA = {
    "data_pds": ["SinglePhoton"],
    "loose": "passGammaCR_Selection && gamma_trigWindowFired",
    "extra_branches": ["passGammaCR_Selection", "HZZ2l2nu_isGammaCR", "gamma_pt", "gamma_eta", "gamma_phi",
                       "gamma_trigWindowFired", "gamma_trigWindowPath", "gamma_prescaleWeight"],
}

BRANCHES = [
    "run", "luminosityBlock", "event",
    # selection
    "passZZ2l2nuSelection", "passZZ2l2nu_emuCR_Selection", "HZZ2l2nu_isEMuCR",
    "HZZ2l2qNu_isELE", "HZZ2l2qNu_cutOppositeChargeFlag", "HZZ2l2nu_nExtraLep", "HZZ2l2nu_nExtraSoftMu", "HZZ2l2nu_nExtraSoftMuPV",
    "HZZ2l2nu_passDeltaRl1l2", "HZZ2l2nu_ifVBF",
    "Triggers_HZZ2l2nu_SR", "Triggers_HZZ2l2nu_EMuCR", "Triggers_HZZ2l2nu_SingleLep",
    # weights (MC only; dropped for data)
    "overallEventWeight", "genWeight", "puWeight", "L1PreFiringWeight_Nom",
    # kinematics
    "massZ1", "pTZ1", "etaZ1", "phiZ1",
    "pTL1", "etaL1", "phiL1", "pTL2", "etaL2", "phiL2", "DeltaRl1l2",
    "pT_MET", "phi_MET", "PuppiMET_pt", "PuppiMET_phi",
    "HZZ2l2nu_ZZmT", "HZZ2l2nu_ZZpT", "HZZ2l2nu_dPhi_ZMET", "HZZ2l2nu_minDPhi_METAK4",
    "HZZ2l2qNu_nJets", "HZZ2l2qNu_nLooseBtagJets", "HZZ2l2qNu_nMediumBtagJets",
    "HZZ2l2qNu_nTightBtagJets",
    "pTj1", "etaj1", "phij1", "pTj2", "etaj2", "phij2",
    "HZZ2l2nu_VBFdijet_mass", "HZZ2l2nu_VBFdEta_jj", "HZZ2l2nu_VBFIndexJet1", "HZZ2l2nu_VBFIndexJet2",
    "PV_npvsGood",
]
MC_ONLY = {"overallEventWeight", "genWeight", "puWeight", "L1PreFiringWeight_Nom"}

ROOT.gInterpreter.Declare("""
#include <set>
#include <tuple>
std::set<std::tuple<UInt_t, UInt_t, ULong64_t>> hzz_seen;
bool hzz_first(UInt_t r, UInt_t l, ULong64_t e) { return hzz_seen.insert({r, l, e}).second; }
""")


def files_of(skimdir, sample):
    return sorted(glob.glob(os.path.join(skimdir, sample, "*_Skim.root")))


def snapshot(files, out, branches, dedup=False, loose=LOOSE):
    df = ROOT.RDataFrame("Events", files)
    avail = set(str(c) for c in df.GetColumnNames())
    cols = [b for b in branches if b in avail]
    missing = [b for b in branches if b not in avail]
    sel = df.Filter(loose, "loose")
    if dedup:
        sel = sel.Filter("hzz_first(run, luminosityBlock, event)", "dedup")
    n_in = df.Count()
    n_out = sel.Count()
    opts = ROOT.RDF.RSnapshotOptions()
    opts.fCompressionAlgorithm = ROOT.RCompressionSetting.EAlgorithm.kLZ4
    opts.fCompressionLevel = 4
    sel.Snapshot("Events", out + ".tmp", cols, opts)
    os.replace(out + ".tmp", out)
    return n_in.GetValue(), n_out.GetValue(), missing


def sum_genweights(files):
    return float(ROOT.RDataFrame("Runs", files).Sum("genEventSumw").GetValue())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["dilepton", "gamma"], default="dilepton")
    ap.add_argument("--skimdir", required=True, help="skim area: <skimdir>/<sample>/*_Skim.root")
    ap.add_argument("--outdir", required=True, help="where the slim files and sumw.json go")
    ap.add_argument("--data_file", default="data_2017.root", help="name of the merged data file in --outdir")
    ap.add_argument("--only", nargs="*", help="process only these sample names (or 'data')")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    data_pds, loose, branches = DATA_PDS, LOOSE, BRANCHES
    if args.mode == "gamma":
        data_pds, loose, branches = GAMMA["data_pds"], GAMMA["loose"], BRANCHES + GAMMA["extra_branches"]

    os.makedirs(args.outdir, exist_ok=True)
    sumw_file = os.path.join(args.outdir, "sumw.json")
    sumw = json.load(open(sumw_file)) if os.path.exists(sumw_file) else {}

    # "_metadata" (written by slurm_setup.py) is not a sample
    samples = sorted(d for d in os.listdir(args.skimdir)
                     if os.path.isdir(os.path.join(args.skimdir, d)) and not d.startswith("_"))
    mc = [s for s in samples if s not in data_pds]

    # ---- data: one deduplicated file, single-threaded (the dedup set is not thread safe)
    if not args.only or "data" in args.only:
        out = os.path.join(args.outdir, args.data_file)
        if os.path.exists(out) and not args.force:
            print("[skip] %s exists" % args.data_file, flush=True)
        else:
            ROOT.DisableImplicitMT()
            files = []
            for pd in data_pds:
                f = files_of(args.skimdir, pd)
                print("data %-15s %4d files" % (pd, len(f)), flush=True)
                files += f
            n_in, n_out, missing = snapshot(files, out, [b for b in branches if b not in MC_ONLY], dedup=True, loose=loose)
            print("%s: %d input events -> %d loose & unique" % (args.data_file, n_in, n_out), flush=True)
            if missing:
                print("  missing branches:", missing, flush=True)

    # ---- MC: one file per sample, multithreaded
    ROOT.EnableImplicitMT(args.threads)
    for s in mc:
        if args.only and s not in args.only:
            continue
        out = os.path.join(args.outdir, s + ".root")
        if os.path.exists(out) and s in sumw and not args.force:
            print("[skip] %s" % s, flush=True)
            continue
        files = files_of(args.skimdir, s)
        if not files:
            print("[warn] no files for %s" % s, flush=True)
            continue
        sw = sum_genweights(files)
        n_in, n_out, missing = snapshot(files, out, branches, loose=loose)
        sumw[s] = {"genEventSumw": sw, "nfiles": len(files), "events_in": n_in, "events_loose": n_out}
        json.dump(sumw, open(sumw_file, "w"), indent=1, sort_keys=True)
        print("%-75s files=%4d sumw=%.6g in=%d loose=%d" % (s, len(files), sw, n_in, n_out), flush=True)
        if missing:
            print("  missing branches:", missing, flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
