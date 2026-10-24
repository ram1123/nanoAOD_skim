#!/usr/bin/env python3
"""Write the 2017 v15 sample list and control-region cut files for myControlPlots.C,
with the event selection synced to AN-2016/325 (see README.md).

Sample list format (myControlPlots.C loadSamples): 8 whitespace-separated columns
  name  file  xsec[pb]  otherscale  nMCevents  nMCNegEvents  colorcode  stackit
myControlPlots scales MC by xsec * lumi * otherscale / (nMCevents - 2*nMCNegEvents).
Here nMCevents=1, nMCNegEvents=0 and otherscale = 1/sum(genEventSumw), so the scale is
xsec * lumi / sum(genEventSumw) -- the normalisation used in higgs_combine/.
Together with gMCWeight = "overallEventWeight" in the RunMacro this is the correct
genWeight-aware normalisation (an event count would be wrong: genWeight can be ~100).

Usage: python make_config.py [--inputs DIR]   (needs sumw.json from ../scripts/analysis/prepare_inputs.py)
"""
import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = "/depot/cms/users/shar1172/HZZ2l2nu_skim/controlplots_2017_AN325/inputs"
LUMI_PB = 41480.0  # UL2017 golden JSON

# (group, sample dir, xsec [pb], colorcode). Stack order = list order (first = bottom).
# Cross sections: taken from this repo's 2018 sample lists where the process exists there
# (HZZ2l2nu_Samples13TeV_ptz_dy_v15.txt and its commented DY lines); the rest are XSDB / GenXSecAnalyzer values.
MC = [
    ("Triboson", "WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8", 0.1707, 34),
    ("Triboson", "WZZ_TuneCP5_13TeV-amcatnlo-pythia8", 0.05709, 34),
    ("Triboson", "ZZZ_TuneCP5_13TeV-amcatnlo-pythia8", 0.01476, 34),
    ("TTV", "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8", 0.2163, 30),
    ("TTV", "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8", 0.2439, 30),
    ("TTV", "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8", 3.697, 30),
    ("TTV", "TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8", 2.967, 30),
    # W+jets: HT-binned LO samples (positive weights). The inclusive amc@NLO sample has only ~250 selected
    # events with weights of +-145 each, which put single-event spikes / negative holes in the stack.
    # No HT<100 bin exists in the 2017 skim (negligible for genuine-dilepton regions).
    ("W+jets", "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8", 1244.0, 46),
    ("W+jets", "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8", 337.8, 46),
    ("W+jets", "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8", 44.93, 46),
    ("W+jets", "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8", 11.19, 46),
    ("W+jets", "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8", 4.926, 46),
    ("W+jets", "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8", 1.152, 46),
    ("W+jets", "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8", 0.02646, 46),
    ("SingleTop", "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8", 3.549, 4),
    ("SingleTop", "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8", 119.7, 4),
    ("SingleTop", "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8", 67.93, 4),
    ("SingleTop", "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", 35.85, 4),
    ("SingleTop", "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", 35.85, 4),
    ("Diboson", "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8", 11.09, 6),
    ("Diboson", "WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8", 49.997, 6),
    ("Diboson", "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8", 5.257, 6),
    ("Diboson", "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8", 6.419, 6),
    ("Diboson", "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8", 0.9738, 6),
    ("Diboson", "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8", 3.676, 6),
    ("Diboson", "GluGluToContinToZZTo2e2nu_TuneCP5_13TeV-mcfm701-pythia8", 0.01719, 6),
    ("Diboson", "GluGluToContinToZZTo2mu2nu_TuneCP5_13TeV-mcfm701-pythia8", 0.01719, 6),
    ("TTbar", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8", 88.29, 7),
    ("Z+jets", "DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8", 20657.0, 2),
    ("Z+jets", "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8", 6198.0, 2),
]
# Not used: QCD, GJets, ZJetsToNuNu, ZNuNuGJets (no genuine dilepton pairs), ZGToLLG / WGToLNuG
# (overlap with DY / W+jets without photon overlap removal), H->ZZ signal samples.

# ---------------------------------------------------------------- control regions
# AN-2016/325 Sec. 4.4 pre-selection, on top of the skim (passZZ2l2nuSelection = 2 tight same-flavour
# leptons pT > 25, no extra loose lepton / soft muon, AN VBF tag). UL / v15 object IDs are kept.
# No charge requirement (AN 4.3). The b-veto uses the UParTAK4 LOOSE WP (AN: CSVv2 loose).
# The ee / mumu signal region (|m_ll - mZ| < 15, MET > 125) is blinded: no region below covers it.
MZ = 91.1876
ONZ = "(abs(massZ1-%g)<15)" % MZ
AN_PRESEL = [
    ("lepton pt", "((pTL1>25)&&(pTL2>25))"),
    ("pTZ1 > 55", "(pTZ1>55)"),
    # AN 4.4.4 soft-muon veto (soft ID, pT > 3). Stored by the skim, not cut there. On 2017 VBF signal it is
    # ~89% efficient (AN: 96-99%); the PV-compatible variant HZZ2l2nu_nExtraSoftMuPV==0 is ~95%.
    ("no extra soft muon", "(HZZ2l2nu_nExtraSoftMu==0)"),
    ("min dPhi(jet,MET) > 0.5", "(HZZ2l2nu_minDPhi_METAK4>0.5)"),
    ("dPhi(Z,MET) > 0.5", "(HZZ2l2nu_dPhi_ZMET>0.5)"),
]
BVETO = [("No loose b-tag", "(HZZ2l2qNu_nLooseBtagJets==0)")]
BTAG = [(">= 1 loose b-tag", "(HZZ2l2qNu_nLooseBtagJets>=1)")]
LL = [("2l2nu event selection", "(passZZ2l2nuSelection==1)"), ("trigger", "(Triggers_HZZ2l2nu_SR==1)")]
EMU = [("emu selection", "(passZZ2l2nu_emuCR_Selection==1)"), ("trigger", "(Triggers_HZZ2l2nu_EMuCR==1)")]
# alpha-method sidebands (AN 5.1): 40 < m_ll < 70 or 110 < m_ll < 200, MET > 70, >= 1 b-tag
ALPHA_SB = [("m_ll sidebands", "(((massZ1>40)&&(massZ1<70))||((massZ1>110)&&(massZ1<200)))"),
            ("pT_MET > 70", "(pT_MET>70)")]
# region -> {flavour: cuts}; every region is split into the categories in JETCATS_OF[region]
REGIONS = {
    "LowMET": {f: LL + c + AN_PRESEL + BVETO + [("on-Z", ONZ), ("pT_MET < 125", "(pT_MET<125)")]
               for f, c in (("ee", [("ee", "(HZZ2l2qNu_isELE==1)")]), ("mumu", [("mumu", "(HZZ2l2qNu_isELE==0)")]))},
    "mllSideband": {f: LL + c + AN_PRESEL + BVETO + [("off-Z", "(!%s)" % ONZ), ("pT_MET > 100", "(pT_MET>100)")]
                    for f, c in (("ee", [("ee", "(HZZ2l2qNu_isELE==1)")]), ("mumu", [("mumu", "(HZZ2l2qNu_isELE==0)")]))},
    # emu events passing the full AN selection, on-Z: N_emu^in of the alpha method
    "emu": {"": EMU + AN_PRESEL + BVETO + [("on-Z", ONZ), ("pT_MET > 125", "(pT_MET>125)")]},
    "alphaSB": {"ee": LL + [("ee", "(HZZ2l2qNu_isELE==1)")] + AN_PRESEL + BTAG + ALPHA_SB,
                "mumu": LL + [("mumu", "(HZZ2l2qNu_isELE==0)")] + AN_PRESEL + BTAG + ALPHA_SB,
                "emu": EMU + AN_PRESEL + BTAG + ALPHA_SB},
}
JETCATS = {
    "0jet": [("VBF category", "(HZZ2l2nu_ifVBF==0)"), ("Zero jet category", "(HZZ2l2qNu_nJets==0)")],
    "ge1jet": [("VBF category", "(HZZ2l2nu_ifVBF==0)"), ("ge1 jet category", "(HZZ2l2qNu_nJets>=1)")],
    "vbf": [("VBF category", "(HZZ2l2nu_ifVBF==1)")],
    "incl": [],
}
JETCATS_OF = {"LowMET": ["0jet", "ge1jet", "vbf"], "mllSideband": ["0jet", "ge1jet", "vbf"],
              "emu": ["0jet", "ge1jet", "vbf"], "alphaSB": ["incl"]}


# ---------------------------------------------------------------- plotted variables
# (branch, xlabel, slog) ; slog only sets the y-axis label format (1: %.1f, 3: GeV, 6: rad).
# Ranges per region (lo, hi, nbins); jet variables only in categories with jets.
#                 var                         xlabel                             slog  LowMET            mllSideband       emu               alphaSB
VARS = [
    ("massZ1", "m_{ll} (GeV)", 3, (76, 106, 30), (10, 200, 38), (76, 106, 10), (40, 200, 32)),
    ("pTZ1", "p_{T}^{ll} (GeV)", 3, (55, 455, 40), (55, 455, 20), (55, 455, 10), (55, 455, 20)),
    ("etaZ1", "#eta^{ll}", 1, (-3, 3, 30), (-3, 3, 15), (-3, 3, 10), (-3, 3, 15)),
    ("phiZ1", "#phi^{ll}", 6, (-3.15, 3.15, 21), (-3.15, 3.15, 14), (-3.15, 3.15, 7), (-3.15, 3.15, 14)),
    ("pTL1", "Leading lepton p_{T} (GeV)", 3, (25, 325, 30), (25, 325, 15), (25, 325, 10), (25, 325, 15)),
    ("pTL2", "Subleading lepton p_{T} (GeV)", 3, (25, 175, 30), (25, 175, 15), (25, 175, 10), (25, 175, 15)),
    ("etaL1", "Leading lepton #eta", 1, (-2.5, 2.5, 25), (-2.5, 2.5, 10), (-2.5, 2.5, 10), (-2.5, 2.5, 10)),
    ("etaL2", "Subleading lepton #eta", 1, (-2.5, 2.5, 25), (-2.5, 2.5, 10), (-2.5, 2.5, 10), (-2.5, 2.5, 10)),
    ("DeltaRl1l2", "#DeltaR(l_{1}, l_{2})", 1, (0, 3.5, 35), (0, 3.5, 14), (0, 3.5, 7), (0, 3.5, 14)),
    ("pT_MET", "p_{T}^{miss} (GeV)", 3, (0, 125, 25), (100, 500, 20), (125, 525, 8), (70, 470, 20)),
    ("phi_MET", "#phi^{miss}", 6, (-3.15, 3.15, 21), (-3.15, 3.15, 14), (-3.15, 3.15, 7), (-3.15, 3.15, 14)),
    ("HZZ2l2nu_ZZmT", "m_{T}^{ZZ} (GeV)", 3, (170, 470, 30), (100, 700, 24), (200, 800, 12), (150, 750, 24)),
    ("HZZ2l2nu_dPhi_ZMET", "#Delta#phi(Z, p_{T}^{miss})", 6, (0.5, 3.2, 27), (0.5, 3.2, 14), (0.5, 3.2, 9), (0.5, 3.2, 14)),
    ("PV_npvsGood", "Number of good primary vertices", 1, (0, 80, 40), (0, 80, 20), (0, 80, 10), (0, 80, 20)),
]
JETVARS = [  # >= 1 jet, VBF and inclusive categories only
    ("HZZ2l2nu_minDPhi_METAK4", "min #Delta#phi(jet, p_{T}^{miss})", 6, (0.5, 3.2, 27), (0.5, 3.2, 14), (0.5, 3.2, 9), (0.5, 3.2, 14)),
    ("HZZ2l2qNu_nJets", "Number of jets", 1, (0, 10, 10), (0, 10, 10), (0, 10, 10), (0, 10, 10)),
    ("HZZ2l2qNu_nLooseBtagJets", "Number of loose b-tagged jets", 1, (0, 5, 5), (0, 5, 5), (0, 5, 5), (0, 5, 5)),
    ("pTj1", "Leading jet p_{T} (GeV)", 3, (30, 430, 40), (30, 630, 20), (30, 630, 10), (30, 630, 20)),
    ("etaj1", "Leading jet #eta", 1, (-4.7, 4.7, 47), (-4.7, 4.7, 20), (-4.7, 4.7, 10), (-4.7, 4.7, 20)),
]
VBFVARS = [  # VBF category only (AN-2016/325 Figs. 12-14)
    ("HZZ2l2nu_VBFdijet_mass", "m_{jj} (GeV)", 3, (500, 2500, 20), (500, 2500, 10), (500, 2500, 5), (500, 2500, 10)),
    ("HZZ2l2nu_VBFdEta_jj", "|#Delta#eta_{jj}|", 1, (4, 9, 20), (4, 9, 10), (4, 9, 5), (4, 9, 10)),
    ("etaj2", "Subleading jet #eta", 1, (-4.7, 4.7, 24), (-4.7, 4.7, 12), (-4.7, 4.7, 6), (-4.7, 4.7, 12)),
]
REGION_INDEX = {"LowMET": 3, "mllSideband": 4, "emu": 5, "alphaSB": 6}

# ---------------------------------------------------------------- single-photon (gamma+jets) CR
# Skim: --channels 1g (nanoAOD_skim branch HZZ_Analysis_2l2nu_dev_gammaCR); inputs: prepare_inputs.py --mode gamma.
# Data = SinglePhoton, HLT_Photon ladder: each event used only via the path of its photon-pT window
# (gamma_trigWindowFired) and weighted by gamma_prescaleWeight (RunMacro sets gDataWeight).
INPUTS_GAMMA = "/depot/cms/users/shar1172/HZZ2l2nu_skim/controlplots_2017/inputs_gamma"
# Cross sections: XSDB / GenXSecAnalyzer values for the UL17 samples -- PLEASE CHECK.
# QCD: HT-binned only (the EMEnriched pT-binned samples overlap with it). GJets and QCD overlap in
# prompt-photon production (no overlap removal applied), as in the AN's gamma-sample composition.
GAMMA_MC = [
    ("Wgamma", "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8", 412.7, 30),
    ("Zgamma", "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8", 55.48, 38),
    ("Znunugamma", "ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-amcatnloFXFX-pythia8", 2.362, 9),
    ("Znunugamma", "ZNuNuGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8", 0.2828, 9),
    ("TopGamma", "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8", 3.697, 7),
    ("TopGamma", "TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8", 2.967, 7),
    ("W+jets", "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8", 1244.0, 46),
    ("W+jets", "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8", 337.8, 46),
    ("W+jets", "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8", 44.93, 46),
    ("W+jets", "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8", 11.19, 46),
    ("W+jets", "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8", 4.926, 46),
    ("W+jets", "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8", 1.152, 46),
    ("W+jets", "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8", 0.02646, 46),
    ("QCD", "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8", 23590000.0, 28),
    ("QCD", "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8", 1555000.0, 28),
    ("QCD", "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8", 324500.0, 28),
    ("QCD", "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8", 30310.0, 28),
    ("QCD", "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8", 6444.0, 28),
    ("QCD", "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8", 1127.0, 28),
    ("QCD", "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8", 109.8, 28),
    ("QCD", "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8", 21.98, 28),
    ("gamma+jets", "GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8", 18540.0, 5),
    ("gamma+jets", "GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8", 8644.0, 5),
    ("gamma+jets", "GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8", 2183.0, 5),
    ("gamma+jets", "GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8", 260.2, 5),
    ("gamma+jets", "GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8", 86.58, 5),
]
GAMMA_CUTS = [  # same AN-2016/325 jet / b-veto / dPhi cuts as the dilepton regions (photon = Z)
    ("gamma CR selection", "(passGammaCR_Selection==1)"),
    ("trigger of the photon-pT window", "(gamma_trigWindowFired==1)"),
    ("pT_gamma > 55", "(gamma_pt>55)"),
    ("min dPhi(jet,MET) > 0.5", "(HZZ2l2nu_minDPhi_METAK4>0.5)"),
    ("dPhi(gamma,MET) > 0.5", "(HZZ2l2nu_dPhi_ZMET>0.5)"),
    ("No loose b-tag", "(HZZ2l2qNu_nLooseBtagJets==0)"),
]
# (branch, xlabel, slog, (lo, hi, nbins))
GAMMA_VARS = [
    ("gamma_pt", "Photon p_{T} (GeV)", 3, (55, 1055, 50)),
    ("gamma_eta", "Photon #eta", 1, (-1.5, 1.5, 30)),
    ("gamma_phi", "Photon #phi", 6, (-3.15, 3.15, 21)),
    ("pT_MET", "p_{T}^{miss} (GeV)", 3, (0, 500, 50)),
    ("phi_MET", "#phi^{miss}", 6, (-3.15, 3.15, 21)),
    ("HZZ2l2nu_ZZmT", "m_{T} (#gamma as Z, GeV)", 3, (150, 1150, 50)),
    ("HZZ2l2nu_dPhi_ZMET", "#Delta#phi(#gamma, p_{T}^{miss})", 6, (0, 3.2, 32)),
    ("PV_npvsGood", "Number of good primary vertices", 1, (0, 80, 40)),
]
GAMMA_JETVARS = [
    ("HZZ2l2nu_minDPhi_METAK4", "min #Delta#phi(jet, p_{T}^{miss})", 6, (0, 3.2, 32)),
    ("HZZ2l2qNu_nJets", "Number of jets", 1, (0, 10, 10)),
    ("HZZ2l2qNu_nLooseBtagJets", "Number of loose b-tagged jets", 1, (0, 5, 5)),
    ("pTj1", "Leading jet p_{T} (GeV)", 3, (30, 630, 30)),
    ("etaj1", "Leading jet #eta", 1, (-4.7, 4.7, 47)),
]


def write_vars(path):
    """controlplotvars_2017.h with one plotVar_t array per region x {0jet, jets}."""
    lines = ["#ifndef _CONTROLPLOTVARS_2017_H", "#define _CONTROLPLOTVARS_2017_H", "",
             '#include "plotvar_t.h"', "",
             "// Generated by make_config.py -- edit VARS / JETVARS there, not here.",
             "// Array name: controlplotvars_2017_<region>_<0jet|jets>; RunMacro_2017_v15.C picks it.", ""]
    for region, idx in REGION_INDEX.items():
        for kind, vlist in (("0jet", VARS), ("jets", VARS + JETVARS), ("vbf", VARS + JETVARS + VBFVARS)):
            lines.append("const plotVar_t controlplotvars_2017_%s_%s[] = {" % (region, kind))
            for v in vlist:
                lo, hi, nb = v[idx]
                lines.append('    {"%s", %g, %g, %d, %d, "%s", "%s", %g, %g, %d, 1},' % (
                    v[0], lo, hi, nb, v[2], v[1], v[0], lo, hi, nb))
            lines.append('    {"", 0.0, 0.0, 0, 0, "", "", 0., 0., 0, 0}};')
            lines.append("")
    for kind, vlist in (("0jet", GAMMA_VARS), ("jets", GAMMA_VARS + GAMMA_JETVARS)):
        lines.append("const plotVar_t controlplotvars_2017_gammaCR_%s[] = {" % kind)
        for v in vlist:
            lo, hi, nb = v[3]
            lines.append('    {"%s", %g, %g, %d, %d, "%s", "%s", %g, %g, %d, 1},' % (
                v[0], lo, hi, nb, v[2], v[1], v[0], lo, hi, nb))
        lines.append('    {"", 0.0, 0.0, 0, 0, "", "", 0., 0., 0, 0}};')
        lines.append("")
    lines.append("#endif // _CONTROLPLOTVARS_2017_H")
    open(path, "w").write("\n".join(lines) + "\n")


def write_cuts(path, cuts):
    with open(path, "w") as f:
        f.write("#  Name                     Cut string   (generated by make_config.py)\n")
        f.write("##########################################\n")
        for name, cut in cuts:
            f.write('"%s"\t"%s"\n' % (name, cut))


def write_samples(out, inputs, mclist):
    sumw = json.load(open(os.path.join(inputs, "sumw.json")))
    with open(out, "w") as f:
        f.write("# name  file  xsec[pb]/lumi[pb-1]  otherscale(=1/sumGenWeight)  nMCevents  nMCNegEvents  colorcode  stackit\n")
        f.write("# generated by make_config.py; MC scale = xsec * lumi / sum(genEventSumw), event weight = overallEventWeight\n")
        f.write("data %s %.1f 1 1 0 1 0\n" % (os.path.join(inputs, "data_2017.root"), LUMI_PB))
        for group, s, xs, col in mclist:
            if s not in sumw:
                print("WARNING: %s missing from %s/sumw.json -- not written" % (s, inputs))
                continue
            f.write("%s %s %g %.10g 1 0 %d 1\n" % (group, os.path.join(inputs, s + ".root"), xs,
                                                   1.0 / sumw[s]["genEventSumw"], col))
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default=INPUTS)
    ap.add_argument("--inputs_gamma", default=INPUTS_GAMMA)
    args = ap.parse_args()

    os.makedirs(os.path.join(HERE, "cuts"), exist_ok=True)
    ncut = 0
    for region, flavs in REGIONS.items():
        for flav, cuts in flavs.items():
            for jc in JETCATS_OF[region]:
                tag = "_".join(x for x in (region, jc, flav) if x)
                write_cuts(os.path.join(HERE, "cuts", "HZZ2l2nu_2017_Cuts_%s.txt" % tag), cuts + JETCATS[jc])
                ncut += 1
    for jc in ("0jet", "ge1jet", "vbf"):
        jcuts = JETCATS[jc]
        write_cuts(os.path.join(HERE, "cuts", "HZZ2l2nu_2017_Cuts_gammaCR_%s.txt" % jc), GAMMA_CUTS + jcuts)
        ncut += 1
    print("wrote %d cut files in cuts/" % ncut)

    write_vars(os.path.join(HERE, "controlplotvars_2017.h"))
    print("wrote controlplotvars_2017.h")

    write_samples(os.path.join(HERE, "HZZ2l2nu_Samples13TeV_2017_v15.txt"), args.inputs, MC)
    if os.path.exists(os.path.join(args.inputs_gamma, "sumw.json")):
        write_samples(os.path.join(HERE, "HZZ2l2nu_Samples13TeV_2017_v15_gamma.txt"), args.inputs_gamma, GAMMA_MC)
    else:
        print("no %s/sumw.json yet (run prepare_inputs.py --mode gamma) -- gamma sample list not written" % args.inputs_gamma)


if __name__ == "__main__":
    main()
