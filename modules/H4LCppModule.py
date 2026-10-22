from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.modules.common.met_phi_correction import METPhiCorrector, Campaign
import ROOT
import yaml
import json
import os
import re
from collections import OrderedDict

from modules.Helper import *
from modules.METFilters import passFilters

ROOT.PyConfig.IgnoreCommandLineOptions = True


class HZZAnalysisCppProducer(Module):

    #def __init__(self, year, cfgFile, isMC, isFSR, cutFlowJSONFile, channels, DEBUG=False, corrector=None):
    def __init__(self, year, cfgFile, isMC, isFSR, cutFlowJSONFile, channels, DEBUG=False,
                 preselection=None, maxEntries=None, era=None):
        self.loadLibraries()
        self.year = year
        self.era = era  # "2016preVFP" / "2016postVFP"; None for other years
        self.isMC = isMC
        self.cutFlowJSONFile = cutFlowJSONFile
        self.channels = channels # choices=["all", "4l", "2l2q", "2l2v", "1g"],
        self.DEBUG = DEBUG
        self.cfgFile = cfgFile
        self.cfg = self._load_config(cfgFile)
        # 1g = single-photon (gamma+jets) CR; needs a GammaCR block in the year's YAML
        self.doGammaCR = self.channels in ("1g", "all") and "GammaCR" in self.cfg
        if self.channels == "1g" and not self.doGammaCR:
            raise KeyError("--channels 1g needs a GammaCR block in {}".format(cfgFile))
        # Event preselection applied by PostProcessor before this module (None = none),
        # and the -n limit; both are only used to label / count the cut flow.
        self.preselection = preselection
        self.maxEntries = maxEntries if maxEntries else None
        self.nInputEntries = 0

        # Trigger lists resolved once to branch names: evaluating every "event.HLT_x"
        # string with eval() on every event was ~40% of analyze(). Each distinct path
        # is read once per event and the per-list OR is built from those values.
        self.triggerLists = OrderedDict()
        for channel in self.cfg['TriggerChannels']:
            paths = []
            for expr in (self.cfg[channel] or []):
                m = re.match(r"^\s*event\.(\w+)\s*$", expr)
                if not m:
                    raise ValueError("Trigger entry '{}' in {} ({}) must be of the form event.HLT_xxx".format(
                        expr, channel, cfgFile))
                paths.append(m.group(1))
            self.triggerLists[channel] = paths
        self.triggerPaths = sorted({p for paths in self.triggerLists.values() for p in paths})

        # gamma CR: exclusive photon-pT windows, each served by one (mostly prescaled) path:
        # [(ptMin, path, data weight = FullLumi / EffLumi)], ascending in ptMin
        self.gammaWindows = []
        if self.doGammaCR:
            g = self.cfg['GammaCR']
            self.gammaWindows = sorted((float(ptMin), path, float(g['FullLumi']) / float(effLumi))
                                       for path, ptMin, effLumi in g['TriggerWindows'])
            missing = [w[1] for w in self.gammaWindows if w[1] not in self.triggerPaths]
            if missing:
                raise ValueError("GammaCR.TriggerWindows paths not in Triggers_HZZ2l2nu_GammaCR: {}".format(missing))

        # MET-phi corrector built once per job (its lookup table was rebuilt every event,
        # ~20% of analyze()). No UL campaign for other years -> uncorrected PuppiMET.
        # UL2016: one data campaign (run-dependent), separate MC campaigns for preVFP (APV) / postVFP.
        campaign = {2016: Campaign.UL_2016, 2017: Campaign.UL_2017, 2018: Campaign.UL_2018}.get(int(year))
        if int(year) == 2016 and isMC:
            if era not in ("2016preVFP", "2016postVFP"):
                raise ValueError("2016 MC needs era 2016preVFP or 2016postVFP (got {})".format(era))
            campaign = Campaign.UL_2016_APV if era == "2016preVFP" else Campaign.UL_2016_NONAPV
        self.metPhiCorrector = (METPhiCorrector(campaign=campaign, is_data=not isMC, is_puppi=True)
                                if campaign is not None else None)
        # MELA is only used by the 4l channel; constructing it is very slow and
        # hangs on some batch nodes, so build it only when 4l is actually run.
        doMELA = self.channels in ("all", "4l")
        self.worker = ROOT.H4LTools(self.year, self.DEBUG, doMELA)
        self._initialize_worker(self.cfg)
        self.worker.isFSR = isFSR
        #self.corrector = corrector
        self._initialize_counters()

        # Cut-flow layout, in the order the selection applies the cuts.
        # Section: (JSON key, printed title, channels it is printed for, rows).
        # Row: (worker counter, printed label, kind) with kind
        #   "start" - first row of a chain (no percentage)
        #   "cut"   - selection step, printed with % of the previous cut/start row
        #   "info"  - count that removes no events (flag, sub-count, MET > 100, ...)
        # Each counter is an int attribute on the C++ worker (declared and reset in
        # H4LTools.h, incremented in H4LTools.cc).
        _z1_rows = [
            ("HZZ2l2qNu_cut2l",           "OS same-flavour Z candidate",    "start"),
            ("cut_2mu_cutOppositeCharge", "  of which mumu",                "info"),
            ("HZZ2l2qNu_cutpTl1l2",       "lepton pT",                      "cut"),
            ("HZZ2l2qNu_cutETAl1l2",      "lepton |eta|",                   "cut"),
            ("HZZ2l2qNu_cutmZ1Window",    "|m_ll - mZ| window",             "cut"),
            ("HZZ2l2qNu_cutZ1Pt",         "pT(ll)",                         "cut"),
        ]
        self.cutFlowSections = [
            ("2l2nu_Channel", "2l2nu signal region (ee / mumu)", ("2l2v", "all"), _z1_rows + [
                ("HZZ2l2qNu_cutDeltaRl1l2",   "  DeltaR(l1,l2) > 0.3 (flag)",   "info"),
                ("HZZ2l2nu_cutExtraLepVeto",  "extra-lepton veto",              "cut"),
                ("HZZ2l2nu_cutSoftMuVeto",    "soft-muon veto",                 "cut"),
                ("HZZ2l2nu_cutdPhiJetMET",    "min dPhi(jet, MET) = selected",  "cut"),
                ("HZZ2l2nu_cutMETgT100",      "  MET > 100 GeV",                "info"),
            ]),
            ("2l2nu_Channel_emu_CR", "2l2nu emu control region", ("2l2v", "all"), [
                ("HZZemuCR_cut2l",            "exactly 1 tight e + 1 tight mu", "start"),
                ("HZZemuCR_cutpTl1l2",        "lepton pT",                      "cut"),
                ("HZZemuCR_cutETAl1l2",       "lepton |eta|",                   "cut"),
                ("HZZemuCR_cutmZ1Window",     "|m_emu - mZ| window",            "cut"),
                ("HZZemuCR_cutZ1Pt",          "pT(emu)",                        "cut"),
                ("HZZemuCR_cutDeltaRl1l2",    "  DeltaR(e,mu) > 0.3 (flag)",    "info"),
                ("HZZemuCR_cutExtraLepVeto",  "extra-lepton veto",              "cut"),
                ("HZZemuCR_cutSoftMuVeto",    "soft-muon veto",                 "cut"),
                ("HZZemuCR_cutdPhiJetMET",    "min dPhi(jet, MET) = selected",  "cut"),
                ("HZZemuCR_cutMETgT100",      "  MET > 100 GeV",                "info"),
            ]),
            ("1g_Channel_gamma_CR", "single-photon (gamma+jets) control region", ("1g", "all"), [
                ("HZZgammaCR_cut1g",          "exactly 1 tight barrel photon",  "start"),
                ("HZZgammaCR_cutExtraLepVeto", "no loose lepton",               "cut"),
                ("HZZgammaCR_cutSoftMuVeto",  "no soft muon",                   "cut"),
                ("HZZgammaCR_cutdPhiJetMET",  "min dPhi(jet, MET) = selected",  "cut"),
                ("HZZgammaCR_cutMETgT100",    "  MET > 100 GeV",                "info"),
            ]),
            ("2l2nu_MuonObjects", "Muons: events with >= 2 muons passing (cumulative)", ("2l2v", "all"), [
                ("cut_mu_pt",                 "pT",                             "start"),
                ("cut_mu_eta",                "+ |eta|",                        "cut"),
                ("cut_mu_mediumid",           "+ mediumId",                     "cut"),
                ("cut_mu_isglobal_istracker", "+ global or tracker",            "cut"),
                ("cut_mu_iso",                "+ dxy/dz + isolation",           "cut"),
            ]),
            ("2l2q_Channel", "2l2q", ("2l2q", "all"), _z1_rows + [
                ("cut2l1J",                   "boosted: >= 1 Z->qq fat jet",    "info"),
                ("cut2l2j",                   "resolved: >= 2 jets",            "info"),
                ("cut2l1Jor2j",               "boosted or resolved",            "info"),
            ]),
            ("4l_Channel", "4l", ("4l", "all"),
                [(c + f, c + f, "start" if c == "cut" else "cut")
                 for f in ("4e", "4mu", "2e2mu")
                 for c in ("cut", "cutghost", "cutLepPt", "cutQCD", "cutZZ", "cutm4l")]),
        ]
        # "General" rows: (key, printed label, kind, channels it is printed for)
        self.cutFlowGeneral = [
            ("Total",                       "Total (input entries)",       "start", None),
            ("PassPresel",                  ("Pass preselection" if self.isMC else "Pass golden JSON + preselection")
                                            if self.preselection else
                                            ("Analysed" if self.isMC else "Pass golden JSON"),
                                                                           "cut",   None),
            ("PassTrig",                    "Pass trigger",                "cut",   None),
            ("PassMETFilters",              "Pass MET filters",            "cut",   None),
            ("PassZZ2l2nuSelection",        "Selected: 2l2nu SR",          "info",  ("2l2v", "all")),
            ("PassZZ2l2nu_emuCR_Selection", "Selected: 2l2nu emu CR",      "info",  ("2l2v", "all")),
            ("PassGammaCR_Selection",       "Selected: gamma CR",          "info",  ("1g", "all")),
            ("PassZZ2l2qSelection",         "Selected: 2l2q",              "info",  ("2l2q", "all")),
            ("PassZZSelection",             "Selected: 4l",                "info",  ("4l", "all")),
        ]

        # Histogram: one bin per General row, then one per section row (bin label = key)
        self.cutFlowBinLabels = [g[0] for g in self.cutFlowGeneral] + \
                                [row[0] for sec in self.cutFlowSections for row in sec[3]]
        hist_nBins = len(self.cutFlowBinLabels)
        self.CutFlowTable = ROOT.TH1F('cutFlow', 'cutFlow', hist_nBins, 0, hist_nBins)
        for idx, label in enumerate(self.cutFlowBinLabels):
            self.CutFlowTable.GetXaxis().SetBinLabel(idx + 1, label)


    def loadLibraries(self):
        base_path = os.getenv('CMSSW_BASE') + '/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim'
        yaml_cpp_path = os.path.join(base_path, "external/yaml-cpp")

        # Adding yaml-cpp headers to the include path
        ROOT.gSystem.AddIncludePath("-I%s/include" % yaml_cpp_path)
        libraries = [
            'libmcfm_710.so',
            'libJHUGenMELAMELA.so',
            'libjhugenmela.so',
            'libcollier.so',
        ]
        for lib in libraries:
            fullPath = os.path.join(base_path, 'external/JHUGenMELA/MELA/data/el9_amd64_gcc12', lib)
            ROOT.gSystem.Load(fullPath)

        # Load the yaml-cpp library
        yaml_cpp_lib_path = os.path.join(yaml_cpp_path, "build")
        ROOT.gSystem.Load(os.path.join(yaml_cpp_lib_path, "libyaml-cpp.so"))

        # Load the C++ module
        if "/H4LTools_cc.so" not in ROOT.gSystem.GetLibraries():
            print("Load C++ module")
            if base_path:
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/H4LTools.cc+O" % base_path)
            else:
                base_path = "$CMSSW_BASE//src/PhysicsTools/NanoAODTools"
                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gROOT.ProcessLine(
                    ".L %s/include/H4LTools.h" % base_path)

    def _load_config(self, cfgFile):
        with open(cfgFile, 'r') as ymlfile:
            return yaml.safe_load(ymlfile)

    def _initialize_worker(self, cfg):
        self.worker.InitializeElecut(*self._get_nested_values(cfg['Electron'], [
            'pTcut', 'Etacut', 'Loosedxycut', 'Loosedzcut',
            'Isocut', ['BDTWP', 'LowEta', 'LowPT'], ['BDTWP', 'MedEta', 'LowPT'],
            ['BDTWP', 'HighEta', 'LowPT'], ['BDTWP', 'LowEta', 'HighPT'],
            ['BDTWP', 'MedEta', 'HighPT'], ['BDTWP', 'HighEta', 'HighPT']
        ]))

        self.worker.InitializeEleMuCleaning(*self._get_nested_values(cfg['Electron'], [
            ['MuonCleaning', 'dR'], ['MuonCleaning', 'MuPt']
        ]))

        self.worker.InitializeMucut(*self._get_nested_values(cfg['Muon'], [
            'pTcut', 'Etacut', 'Loosedxycut', 'Loosedzcut', 'Isocut',
            'Tightdxycut', 'Tightdzcut', 'TightTrackerLayercut', 'TightpTErrorcut',
            'HighPtBound'
        ]))

        self.worker.InitializeFsrPhotonCut(*self._get_nested_values(cfg['FsrPhoton'], [
            'pTcut', 'Etacut', 'Isocut', 'dRlcut', 'dRlOverPtcut'
        ]))

        self.worker.InitializeJetcut(*self._get_nested_values(cfg['Jet'], ['pTcut', 'Etacut', ['UParTAK4_btag', 'Etacut'], 'LepCleanDR']))

        self.worker.InitializeEvtCut(*self._get_nested_values(cfg, ['MZ1cut', 'MZZcut',
                                                                    ['Higgscut', 'down'], ['Higgscut', 'up'],
                                                                    'Zmass', ['MZcut', 'down'], ['MZcut', 'up'],
                                                                    ['Jet','UParTAK4_btag','Loose'], ['Jet','UParTAK4_btag','Medium'], ['Jet','UParTAK4_btag','Tight']]))

        if self.channels == "all" or self.channels == "2l2q":
            self.worker.InitializeHZZ2l2qCut(*self._get_nested_values(cfg['HZZ2l2q'],
                                                                  ['Leading_Lep_pT', 'SubLeading_Lep_pT', 'Lep_eta',
                                                                   ['MZLepcut', 'down'], ['MZLepcut', 'up']]))

        if self.channels in ("all", "2l2v", "1g"):   # 1g reuses the 2l2nu jet / veto / dPhi settings
            self.worker.InitializeHZZ2l2nuCut(*self._get_nested_values(cfg['HZZ2l2nu'],
                                                                   ['Leading_Lep_pT', 'SubLeading_Lep_pT', 'Lep_eta', 'Pt_ll',
                                                                    'M_ll_Window', 'dPhi_jetMET', ['MZLepcut', 'down'], ['MZLepcut', 'up']]))
            self.worker.InitializeVetoLepCut(*self._get_nested_values(cfg['HZZ2l2nu'],
                                                                   [['VetoLep', 'Ele_pT'], ['VetoLep', 'Ele_cutBased'], ['VetoLep', 'Ele_Iso'],
                                                                    ['VetoLep', 'Mu_pT'], ['VetoLep', 'Mu_Iso']]))
            self.worker.InitializeEmuCRCut(*self._get_nested_values(cfg['HZZ2l2nu'],
                                                                   [['emuCR', 'M_ll_Window'], ['emuCR', 'Pt_ll']]))
            # AN-2016/325: charge requirement, soft-muon veto, VBF tag (argument order = H4LTools.h)
            self.worker.InitializeHZZ2l2nuAN(*self._get_nested_values(cfg['HZZ2l2nu'],
                                                                   ['RequireOS', ['VetoLep', 'SoftMu_pT'],
                                                                    ['VBF', 'dEta_jj'], ['VBF', 'M_jj'],
                                                                    ['VBF', 'ZCentrality'], ['VBF', 'CentralJetVeto']]))
        if self.doGammaCR:
            self.worker.InitializeGammaCRCut(*self._get_nested_values(cfg['GammaCR'],
                                                                      ['pTcut', 'IDcut', 'JetCleanDR']))

    def _get_nested_values(self, dictionary, keys):
        # Raise a clear error for a missing key instead of passing a placeholder
        # string into a float Initialize* argument (which fails obscurely later).
        values = []
        for key in keys:
            if isinstance(key, list):
                sub_dict = dictionary
                for sub_key in key:
                    if not isinstance(sub_dict, dict) or sub_key not in sub_dict:
                        raise KeyError("Missing config key: {} (in {})".format(
                            " -> ".join(str(k) for k in key), self.cfgFile))
                    sub_dict = sub_dict[sub_key]
                values.append(sub_dict)
            else:
                if key not in dictionary:
                    raise KeyError("Missing config key: {} (in {})".format(
                        key, self.cfgFile))
                values.append(dictionary[key])
        return values

    def _get_weight_branch(self, event, name, default, warn=None):
        # nanoAOD-tools raises RuntimeError("Unknown branch ...") for an absent
        # branch, which getattr(event, name, default) does NOT catch. Use this so a
        # missing weight producer (e.g. no PU weight for 2022, or a broken re-nano
        # without L1PreFiringWeight_*) degrades to `default` with a one-time warning
        # instead of crashing the job.
        try:
            return getattr(event, name)
        except Exception:
            if warn is not None:
                flag = "_warned_missing_" + name
                if not getattr(self, flag, False):
                    print("WARNING: branch '{}' absent for year {} MC - {} NOT "
                          "applied (using {})".format(name, self.year, warn, default))
                    setattr(self, flag, True)
            return default

    def _initialize_counters(self):
        self.passAllEvts = 0
        self.passtrigEvts = 0
        self.passMETFilters = 0
        self.passZZ4lEvts = 0
        self.passZZ2l2qEvts = 0
        self.passZZ2l2nuEvts = 0
        self.passZZ2l2nu_emuCR_Evts = 0
        self.passGammaCR_Evts = 0

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree)  # initReaders must be called in beginFile
        # Cut-flow "Total": all entries of the input file (within -n), counted before the
        # golden JSON and the preselection that PostProcessor applies ahead of analyze()
        nFile = inputFile.Get("Events").GetEntries()
        self.nInputEntries += min(nFile, self.maxEntries) if self.maxEntries else nFile
        self.out = wrappedOutputTree

        # NanoAOD version from the jet branches: v9 stores Jet_puId (AK4 CHS), v15 does not (AK4 PUPPI)
        self.isNanoV9 = bool(inputTree.GetBranch("Jet_puId"))
        self.worker.nanoAODVersion = 9 if self.isNanoV9 else 15
        print("INFO: NanoAOD v{} jets in {}: {}".format(self.worker.nanoAODVersion, inputFile.GetName(),
              "stored jetId + PU jet ID" if self.isNanoV9 else "recomputed jet ID, no PU jet ID"))

        # A path missing from the HLT menu of this file's era (2016 emu paths) has no branch: it
        # counts as not fired instead of crashing getattr() in analyze().
        self.presentTriggerPaths = [p for p in self.triggerPaths if inputTree.GetBranch(p)]
        missing = sorted(set(self.triggerPaths) - set(self.presentTriggerPaths))
        if missing:
            print("WARNING: trigger paths absent in {} (treated as not fired): {}".format(
                inputFile.GetName(), ", ".join(missing)))

        # Boolean branches for Trigger channels
        for TriggerChannel in self.cfg['TriggerChannels']:
            self.out.branch(TriggerChannel, "O")
        # boolean branches for 4l, 2l2q, 2l2nu channels
        self.out.branch("passZZ4lSelection",  "O")
        self.out.branch("passZZ2l2qSelection",  "O")
        self.out.branch("passZZ2l2nuSelection",  "O")
        self.out.branch("passZZ2l2nu_emuCR_Selection", "O")
        self.out.branch("isBoosted2l2q",  "O")
        self.out.branch("HZZ2l2nu_ifVBF", "O")
        self.out.branch("HZZ2l2qNu_isELE", "O")
        self.out.branch("HZZ2l2qNu_cutOppositeChargeFlag", "O")
        self.out.branch("HZZ2l2nu_isEMuCR", "O")
        # single-photon CR (--channels 1g / all): photon kinematics; pTZ1/etaZ1/phiZ1 hold the
        # photon too (Z1 = photon with m = mZ) so the 2l2nu variables (ZZmT, dPhi_ZMET, jets) apply
        self.out.branch("passGammaCR_Selection", "O")
        self.out.branch("HZZ2l2nu_isGammaCR", "O")
        self.out.branch("gamma_pt", "F")
        self.out.branch("gamma_eta", "F")
        self.out.branch("gamma_phi", "F")
        self.out.branch("gamma_trigWindowFired", "O")   # path of the photon's pT window fired (require it)
        self.out.branch("gamma_trigWindowPath", "I")    # index of that path in GammaCR.TriggerWindows (-1 = none)
        self.out.branch("gamma_prescaleWeight", "F")    # data: FullLumi / EffLumi of that path; MC: 1

        # Branches for leptons related varialbes: 4l, 2l2q, 2l2nu
        self.out.branch("massL1",  "F")
        self.out.branch("pTL1",  "F")
        self.out.branch("etaL1",  "F")
        self.out.branch("phiL1",  "F")
        self.out.branch("massL2",  "F")
        self.out.branch("pTL2",  "F")
        self.out.branch("etaL2",  "F")
        self.out.branch("phiL2",  "F")
        self.out.branch("massL3",  "F")
        self.out.branch("pTL3",  "F")
        self.out.branch("etaL3",  "F")
        self.out.branch("phiL3",  "F")
        self.out.branch("massL4",  "F")
        self.out.branch("pTL4",  "F")
        self.out.branch("etaL4",  "F")
        self.out.branch("phiL4",  "F")
        self.out.branch("DeltaRl1l2",  "F")
        self.out.branch("HZZ2l2nu_passDeltaRl1l2", "O")
        # Branches for 4l channel: ZZ kinematics
        self.out.branch("mass4l",  "F")
        self.out.branch("pT4l",  "F")
        self.out.branch("eta4l",  "F")
        self.out.branch("phi4l",  "F")

        # Branches for 4l channel: MELA discriminants
        self.out.branch("D_CP",  "F")
        self.out.branch("D_0m",  "F")
        self.out.branch("D_0hp",  "F")
        self.out.branch("D_int",  "F")
        self.out.branch("D_L1",  "F")
        self.out.branch("D_L1Zg",  "F")

        # Branches for 4l channel only; These contains the kinematics of the AK4 jets
        self.out.branch("mj1",  "F")
        self.out.branch("pTj1",  "F")
        self.out.branch("etaj1",  "F")
        self.out.branch("phij1",  "F")

        self.out.branch("mj2",  "F")
        self.out.branch("pTj2",  "F")
        self.out.branch("etaj2",  "F")
        self.out.branch("phij2",  "F")

        # common branches for 4l, 2l2q, 2l2nu channels
        self.out.branch("massZ1",  "F")
        self.out.branch("pTZ1",  "F")
        self.out.branch("etaZ1",  "F")
        self.out.branch("phiZ1",  "F")

        self.out.branch("massZ2",  "F")
        self.out.branch("pTZ2",  "F")
        self.out.branch("etaZ2",  "F")
        self.out.branch("phiZ2",  "F")

        self.out.branch("pT_MET",  "F")
        self.out.branch("phi_MET",  "F")

        # --- MC event weights (all = 1.0 for data) -------------------------------
        # overallEventWeight        = genWeight * puWeight     * L1PreFiringWeight_Nom
        # overallEventWeight_puUp   = genWeight * puWeightUp   * L1PreFiringWeight_Nom
        # overallEventWeight_puDown = genWeight * puWeightDown * L1PreFiringWeight_Nom
        # overallEventWeight_prefire{Up,Down} = genWeight * puWeight * L1PreFiringWeight_{Up,Dn}
        # Each variation shifts exactly ONE factor. Only weight-based systematics
        # that do NOT change the event selection are here (PU, L1 prefiring);
        # JES/JER/unclustered-MET/lepton-scale need a per-variation re-run of the
        # C++ selection and are NOT included.
        #
        # DOWNSTREAM CONTRACT (higgs_combine/make_shapes*.cpp must follow this):
        #   yield = sum(overallEventWeight) * xsec * lumi / genEventSumw
        #   - genEventSumw is the sum of genWeight over ALL generated events, from
        #     the merged Runs tree (preserved by haddnano.py). Do NOT normalise by
        #     the event count or by (nPos - nNeg): those differ from genEventSumw
        #     for samples whose genWeight is not +/-1 (POWHEG MiNLO, aMC@NLO).
        #   - Do NOT multiply by genWeight or puWeight again: they are folded in
        #     here. The standalone genWeight / puWeight / L1PreFiringWeight_* branches
        #     are kept for provenance/debugging only.
        #   - LHEScaleWeight / LHEPdfWeight / PSWeight are deliberately NOT folded
        #     in; apply them downstream as ratios renormalised to LHEScaleSumw /
        #     LHEPdfSumw (scale/PDF) or genEventSumw (PS).
        self.out.branch("overallEventWeight",            "F")
        self.out.branch("overallEventWeight_puUp",        "F")
        self.out.branch("overallEventWeight_puDown",      "F")
        self.out.branch("overallEventWeight_prefireUp",   "F")
        self.out.branch("overallEventWeight_prefireDown", "F")

        # Branches for 2l2q channel
        self.out.branch("massZ2_2j",  "F")
        self.out.branch("phiZ2_2j",  "F")
        self.out.branch("etaZ2_2j",  "F")
        self.out.branch("pTZ2_2j",  "F")
        self.out.branch("EneZ2_2j",  "F")
        self.out.branch("phiZ2_met",  "F")
        self.out.branch("pTZ2_met",  "F")
        self.out.branch("EneZ2_met",  "F")

        # Branches for 2l2nu channel: ZZ kinematics
        self.out.branch("HZZ2l2nu_ZZmT",  "F")
        self.out.branch("HZZ2l2nu_ZZpT",  "F")
        #self.out.branch("Pz_neutrino", "F")

        # Branches for 2l2nu channel: VBF jets and dijet kinematics
        self.out.branch("HZZ2l2qNu_nJets", "I")
        self.out.branch("HZZ2l2nu_nExtraLep", "I")
        self.out.branch("HZZ2l2nu_nExtraSoftMu", "I")    # AN-2016/325 soft-muon veto (pT > 3, soft ID); == 0 required in the skim
        self.out.branch("HZZ2l2nu_nExtraSoftMuPV", "I")  # same, |dz| < 0.5, |dxy| < 0.2 cm

        self.out.branch("HZZ2l2qNu_nTightBtagJets", "I")
        self.out.branch("HZZ2l2qNu_nMediumBtagJets", "I")
        self.out.branch("HZZ2l2qNu_nLooseBtagJets", "I")

        self.out.branch("HZZ2l2nu_VBFIndexJet1",  "I")
        self.out.branch("HZZ2l2nu_VBFIndexJet2",  "I")
        self.out.branch("HZZ2l2nu_minDPhi_METAK4",  "F")
        self.out.branch("HZZ2l2nu_dPhi_ZMET",  "F")  # |dphi(Z, MET)|; stored only, no cut applied

        self.out.branch("HZZ2l2nu_VBFjet1_pT",  "F")
        self.out.branch("HZZ2l2nu_VBFjet1_eta",  "F")
        self.out.branch("HZZ2l2nu_VBFjet1_phi",  "F")
        self.out.branch("HZZ2l2nu_VBFjet1_mass",  "F")

        self.out.branch("HZZ2l2nu_VBFjet2_pT",  "F")
        self.out.branch("HZZ2l2nu_VBFjet2_eta",  "F")
        self.out.branch("HZZ2l2nu_VBFjet2_phi",  "F")
        self.out.branch("HZZ2l2nu_VBFjet2_mass",  "F")

        self.out.branch("HZZ2l2nu_VBFdijet_mass",  "F")
        self.out.branch("HZZ2l2nu_VBFdijet_pT",  "F")
        self.out.branch("HZZ2l2nu_VBFdijet_E",  "F")
        self.out.branch("HZZ2l2nu_VBFdEta_jj",  "F")
        self.out.branch("HZZ2l2nu_VBFdPhi_jj",  "F")
        self.out.branch("HZZ2l2nu_VBFdR_jj",  "F")

        # Branches for 2l2q channel
        self.out.branch("HZZ2l2q_boostedJet_PNScore", "F")
        self.out.branch("HZZ2l2q_boostedJet_Index", "I")
        self.out.branch("HZZ2l2q_resolvedJet1_Index", "I")
        self.out.branch("HZZ2l2q_resolvedJet2_Index", "I")

        # FSR branches for leptons
        self.out.branch("Electron_Fsr_pt",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_eta",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_phi",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Muon_Fsr_pt",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_eta",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_phi",  "F", lenVar = "nMuon_Fsr")

        with open("SyncLepton2018GGH.txt", 'w') as f:
            f.write("Sync data list:"+"\n")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        generalCounts = {
            "Total": self.nInputEntries,
            "PassPresel": self.passAllEvts,
            "PassTrig": self.passtrigEvts,
            "PassMETFilters": self.passMETFilters,
            "PassZZSelection": self.passZZ4lEvts,
            "PassZZ2l2qSelection": self.passZZ2l2qEvts,
            "PassZZ2l2nuSelection": self.passZZ2l2nuEvts,
            "PassZZ2l2nu_emuCR_Selection": self.passZZ2l2nu_emuCR_Evts,
            "PassGammaCR_Selection": self.passGammaCR_Evts,
        }
        cutFlowData = {"General": generalCounts}
        for key, _, _, rows in self.cutFlowSections:
            cutFlowData[key] = {row[0]: int(getattr(self.worker, row[0])) for row in rows}

        def printRows(rows):
            """rows: (label, count, kind); cut rows get % of the previous start/cut row."""
            previous = None
            for label, count, kind in rows:
                pct = ""
                if kind == "cut" and previous:
                    pct = "{:7.2f}%".format(100.0 * count / previous)
                if kind in ("start", "cut"):
                    previous = count
                print("  {:<36}{:>10}  {}".format(label, count, pct))

        def printed(channels):
            return channels is None or self.channels in channels

        print("\n==================== Cut flow ({} channel{}) ====================".format(
            self.channels, "s" if self.channels == "all" else ""))
        print("% = fraction of the previous selection step; indented rows are counts that remove no events")
        print("Preselection: {}\n".format(self.preselection if self.preselection else "none"))
        print("General")
        printRows([(label, generalCounts[key], kind)
                   for key, label, kind, channels in self.cutFlowGeneral if printed(channels)])
        for key, title, channels, rows in self.cutFlowSections:
            if not printed(channels):
                continue
            print("\n" + title)
            printRows([(label, cutFlowData[key][counter], kind) for counter, label, kind in rows])
        # electron-muon cross-cleaning is disabled in H4LTools::LeptonSelection (docs/TODO.md)
        #print("\nElectrons removed by electron-muon cross-cleaning (dR < {:.2f}): {}".format(
        #    self.worker.eleMuCleaningDR, self.worker.nEleMuOverlapRemoved))
        print("\n=================================================================\n")

        # Histogram bins follow self.cutFlowBinLabels
        flat = dict(generalCounts)
        for key, _, _, _ in self.cutFlowSections:
            flat.update(cutFlowData[key])
        for idx, label in enumerate(self.cutFlowBinLabels):
            self.CutFlowTable.SetBinContent(idx + 1, flat[label])

        # Write the cut flow data to a json file
        with open(self.cutFlowJSONFile, "w") as jsonFile:
            json.dump(cutFlowData, jsonFile, indent=5)

        outputFile.cd()
        self.CutFlowTable.Write()
        pass

    # this function gets the pointers to Value and ArrayReaders and sets
    # them in the C++ worker class
    def initReaders(self, tree):
        # self._ttreereaderversion must be set AFTER all calls to
        # tree.valueReader or tree.arrayReader
        self._ttreereaderversion = tree._ttreereaderversion

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail,
        go to next event)"""
        #if event.run != 317292 or event.luminosityBlock != 95 or event.event != 144614313:
        #if event.Electron_pt.GetSize() > 0:
            #return False
        #if event.nElectron != 1 or event.nMuon != 1:
            #return False
        #print("Event electron_pt =", event.Electron_pt)
        #print("Event Muon_pt =", event.Muon_pt)
        #print("Event MET_pt =", event.MET_pt)
        if self.DEBUG:
            print("======       Inside analyze function     ==========")
        # do this check at every event, as other modules might have read
        # further branches
        # if event._tree._ttreereaderversion > self._ttreereaderversion:
        #    self.initReaders(event._tree)
        # do NOT access other branches in python between the check/call to
        # initReaders and the call to C++ worker code
        self.worker.Initialize()
        isMC = self.isMC
        self.worker.SetObjectNum(event.nElectron,event.nMuon,event.nJet,event.nFsrPhoton)
        if isMC:
            self.worker.SetObjectNumGen(event.nGenPart)

        keepIt = False

        # --- MC event weights (see beginFile for the definition) ---
        if isMC:
            _is_run2 = int(self.year) in (2016, 2017, 2018)
            _gen = self._get_weight_branch(event, "genWeight", 1.0)
            _puN = self._get_weight_branch(event, "puWeight", 1.0,
                                           warn="PU reweighting")           # puWeight_UL20XX
            _puU = self._get_weight_branch(event, "puWeightUp", _puN)
            _puD = self._get_weight_branch(event, "puWeightDown", _puN)
            # L1PreFiringWeight_Nom = combined ECAL (2016-17) x muon (2016-18);
            # for 2018 it is driven by the muon term and is NOT ~1.
            _prN = self._get_weight_branch(event, "L1PreFiringWeight_Nom", 1.0,
                                           warn="L1 prefiring" if _is_run2 else None)
            _prU = self._get_weight_branch(event, "L1PreFiringWeight_Up", _prN)
            _prD = self._get_weight_branch(event, "L1PreFiringWeight_Dn", _prN)
            overallEventWeight             = _gen * _puN * _prN
            overallEventWeight_puUp        = _gen * _puU * _prN
            overallEventWeight_puDown      = _gen * _puD * _prN
            overallEventWeight_prefireUp   = _gen * _puN * _prU
            overallEventWeight_prefireDown = _gen * _puN * _prD
        else:
            overallEventWeight = overallEventWeight_puUp = overallEventWeight_puDown = \
                overallEventWeight_prefireUp = overallEventWeight_prefireDown = 1.0

        passedTrig=False
        passedFullSelection=False
        passedZ4lSelection=False
        passedQCDcut=False
        passedZ1LSelection=False
        passedZ4lZ1LSelection=False
        passedZ4lZXCRSelection=False
        passedZXCRSelection=False
        passedFiducialSelection=False
        nZXCRFailedLeptons=0
        self.passAllEvts += 1
        self.CutFlowTable.Fill(0)
        massZ2_2j = -999.
        phiZ2_2j = -999.
        etaZ2_2j = -999.
        pTZ2_2j = -999.
        EneZ2_2j = -999.
        phiZ2_met = -999.
        pTZ2_met = -999.
        EneZ2_met = -999.
        HZZ2l2nu_minDPhi_METAK4 = 999.0
        HZZ2l2nu_dPhi_ZMET = -999.

        HZZ2l2nu_ZZmT = -999.
        HZZ2l2nu_ZZpT = -999.

        HZZ2l2nu_VBFIndexJet1 = -999
        HZZ2l2nu_VBFIndexJet2 = -999

        HZZ2l2nu_VBFjet1_pT = -999.
        HZZ2l2nu_VBFjet1_eta = -999.
        HZZ2l2nu_VBFjet1_phi = -999.
        HZZ2l2nu_VBFjet1_mass = -999.

        HZZ2l2nu_VBFjet2_pT = -999.
        HZZ2l2nu_VBFjet2_eta = -999.
        HZZ2l2nu_VBFjet2_phi = -999.
        HZZ2l2nu_VBFjet2_mass = -999.

        HZZ2l2nu_VBFdijet_mass = -999.
        HZZ2l2nu_VBFdijet_pT = -999.
        HZZ2l2nu_VBFdijet_E = -999.
        HZZ2l2nu_VBFdEta_jj = -999.
        HZZ2l2nu_VBFdPhi_jj = -999.
        HZZ2l2nu_VBFdR_jj = -999.

        pTL1 = -999.
        etaL1 = -999.
        phiL1 = -999.
        massL1 = -999.
        pTL2 = -999.
        etaL2 = -999.
        phiL2 = -999.
        massL2 = -999.
        pTZ1 = -999.
        etaZ1 = -999.
        phiZ1 = -999.
        massZ1 = -999.
        pTZ2 = -999.
        etaZ2 = -999.
        phiZ2 = -999.
        massZ2 = -999.
        D_CP = -999.
        D_0m = -999.
        D_0hp = -999.
        D_int = -999.
        D_L1 = -999.
        D_L1Zg = -999.
        pTL3 = -999.
        etaL3 = -999.
        phiL3 = -999.
        massL3 = -999.
        pTL4 = -999.
        etaL4 = -999.
        phiL4 = -999.
        massL4 = -999.
        pTj1 = -999.
        etaj1 = -999.
        phij1 = -999.
        mj1 = -999.
        pTj2 = -999.
        etaj2 = -999.
        phij2 = -999.
        mj2 = -999.
        pT4l = -999.
        eta4l = -999.
        phi4l = -999.
        mass4l = -999.
        pT4l = -999.
        eta4l = -999.
        phi4l = -999.
        mass4l = -999.
        DeltaRl1l2 = -999.
        passedTrig = False
        fired = {path: bool(getattr(event, path)) for path in self.presentTriggerPaths}
        TriggerMap = {channel: any(fired.get(path, False) for path in paths)
                      for channel, paths in self.triggerLists.items()}

        # If any of the trigger channel from TriggerMap passes, then the event is kept else return keepIt
        for value in TriggerMap.values():
            if value:
                passedTrig = True
                break
        if not passedTrig:
            return keepIt

        self.passtrigEvts += 1
        self.CutFlowTable.Fill(1)

        if passFilters(event, int(self.year), isMC=self.isMC):
            self.passMETFilters += 1
            self.CutFlowTable.Fill(2)
        else:
            return keepIt

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        if self.doGammaCR:
            for xp in Collection(event, "Photon"):
                # argument order must match H4LTools.h SetPhotons
                self.worker.SetPhotons(xp.pt, xp.eta, xp.phi, xp.cutBased, xp.isScEtaEB, xp.electronVeto, xp.pixelSeed)
        jets = Collection(event, "Jet")
        FatJets = Collection(event, "FatJet")
        puppimet = Object(event, "PuppiMET", None)
        corrector = self.metPhiCorrector

        #MET correction for v15
        if corrector is not None:
            corr_pt, corr_phi = corrector(
            puppimet.pt,
            puppimet.phi,
            npv=event.PV_npvs,
            run=event.run
            )
        else:
            if not getattr(self, "_warned_no_metphi", False):
                print("WARNING: no MET-phi correction for year {} - using uncorrected PuppiMET".format(self.year))
                self._warned_no_metphi = True
            corr_pt, corr_phi = puppimet.pt, puppimet.phi

        # for photon in Photons:
        #     # Keep photons if pT > 55, |eta| < 2.5 and skip the transition region of barrel and endcap
        #     if photon.pt > 55 and abs(photon.eta) < 2.5 and not (1.4442 < abs(photon.eta) < 1.566):
        #         return keepIt

        if isMC:
            genparts = Collection(event, "GenPart")
            GenJets = Collection(event, "GenJet")
            for xg in genparts:
                self.worker.SetGenParts(xg.pt)
            for xm in muons:
                self.worker.SetMuonsGen(xm.genPartIdx)
            for xg in GenJets:
                self.worker.SetGenJets(xg.pt, xg.eta, xg.phi, xg.mass)

        for xe in electrons:
            #self.worker.SetElectrons(xe.pt, xe.eta, xe.phi, xe.mass, xe.dxy,
                                      #xe.dz, xe.mvaFall17V2Iso_WP90, xe.pdgId, xe.pfRelIso03_all) #for v9
            self.worker.SetElectrons(xe.pt, xe.eta, xe.phi, xe.mass, xe.dxy, xe.dz, xe.mvaIso_WP90, xe.pfRelIso03_all, xe.pdgId, xe.deltaEtaSC) # for v15
            self.worker.SetElectronsCutBased(xe.cutBased)   # veto-electron ID (Fall17V2 cut-based)
            if self.DEBUG:
                print("Electrons: pT, eta: {}, {}".format(xe.pt, xe.eta))

        for xm in muons:
            self.worker.SetMuons(xm.corrected_pt, xm.eta, xm.phi, xm.mass, xm.isGlobal, xm.isTracker, xm.mediumId,
                                xm.dxy, xm.dz, xm.ptErr, xm.nTrackerLayers, xm.isPFcand,
                                 xm.pdgId, xm.charge, xm.pfRelIso03_all)
            self.worker.SetMuonsSoftId(xm.softId)
            if self.DEBUG:
                print("Muons: pT, eta: {}, {}".format(xm.corrected_pt, xm.eta))

        for xf in fsrPhotons:
            self.worker.SetFsrPhotons(xf.dROverEt2,xf.eta,xf.phi,xf.pt,xf.relIso03)
        
        # for v9
        #for xj in jets:
            #self.worker.SetJets(xj.pt,xj.eta,xj.phi,xj.mass,xj.jetId, xj.btagDeepFlavB, xj.puId)
        # for v15
        for xj in jets:
            # argument order must match H4LTools.h SetJets: chEmEF, chHEF, neEmEF, neHEF, muEF, nConstituents, neMultiplicity, chMultiplicity
            self.worker.SetJets(xj.pt,xj.eta,xj.phi,xj.mass,xj.btagUParTAK4B, xj.chEmEF, xj.chHEF, xj.neEmEF, xj.neHEF, xj.muEF, xj.nConstituents, xj.neMultiplicity, xj.chMultiplicity)
            if self.isNanoV9:
                self.worker.SetJetIDv9(xj.jetId, xj.puId)

        #for xj in FatJets:
            #self.worker.SetFatJets(xj.pt, xj.eta, xj.phi, xj.msoftdrop, xj.jetId, xj.btagDeepB, xj.particleNet_ZvsQCD)
            #self.worker.SetFatJets(xj.pt, xj.eta, xj.phi, xj.msoftdrop)

        #self.worker.SetPuppiMET(puppimet.pt, puppimet.phi, puppimet.sumEt) # for v15
        self.worker.SetPuppiMET(corr_pt, corr_phi, puppimet.sumEt)
        if self.DEBUG:
            print("***** MET: corr_pt, corr_phi, sumEt: {}, {}, {}".format(corr_pt, corr_phi, puppimet.sumEt))
            print("***** MET not corrected: pt, phi, sumEt: {}, {}, {}".format(puppimet.pt, puppimet.phi, puppimet.sumEt))
            print("***** Event branch MET_pt =", event.PuppiMET_pt)

        self.worker.LeptonSelection()
        foundZZCandidate_4l = False    # for 4l
        passZZ4lSelection = False
        foundZZCandidate_2l2q = False # for 2l2q
        passZZ2l2qSelection = False
        foundZZCandidate_2l2nu = False # for 2l2nu
        passZZ2l2nuSelection = False
        foundZZCandidate_2l2nu_emuCR = False   # for 2l2nu emu control region
        passZZ2l2nu_emuCR_Selection = False    # for 2l2nu emu control region
        HZZ2l2nu_ifVBF = False
        HZZ2l2qNu_isELE = False
        HZZ2l2nu_nExtraLep = -1
        HZZ2l2nu_isEMuCR = False
        HZZ2l2qNu_cutOppositeChargeFlag = False
        isBoosted2l2q = False

        if self.worker.GetZ1_2l2qOR2l2nu() and (self.channels == "all"  or self.channels == "2l2v" or self.channels == "2l2q"):  #commented out for now
            if self.channels == "2l2q" or self.channels == "all":
                foundZZCandidate_2l2q = self.worker.ZZSelection_2l2q()
                isBoosted2l2q = self.worker.isBoosted2l2q
                if self.DEBUG: print("isBoosted2l2q: ", isBoosted2l2q)
            if self.channels == "2l2v" or self.channels == "all":
                foundZZCandidate_2l2nu = self.worker.ZZSelection_2l2nu()  #commented out for now
        if (self.channels == "all"  or self.channels == "4l"):
            foundZZCandidate_4l = self.worker.ZZSelection_4l()
        if self.worker.GetZ1_emuCR() and (self.channels == "all"  or self.channels == "2l2v"):
            foundZZCandidate_2l2nu_emuCR = self.worker.ZZSelection_2l2nu()
        # single-photon CR: own trigger; never together with a dilepton candidate (it vetoes leptons)
        foundGammaCR = False
        passGammaCR_Selection = False
        gamma_pt = gamma_eta = gamma_phi = -999.
        gamma_trigWindowFired = False
        gamma_trigWindowPath = -1
        gamma_prescaleWeight = 0.
        if (self.doGammaCR and TriggerMap["Triggers_HZZ2l2nu_GammaCR"]
                and not (foundZZCandidate_4l or foundZZCandidate_2l2q or foundZZCandidate_2l2nu or foundZZCandidate_2l2nu_emuCR)
                and self.worker.GetZ1_gammaCR()):
            foundGammaCR = self.worker.ZZSelection_2l2nu()

        HZZ2l2q_boostedJet_PNScore = self.worker.boostedJet_PNScore
        HZZ2l2q_boostedJet_Index = self.worker.boostedJet_Index
        HZZ2l2q_resolvedJet1_Index = self.worker.resolvedJet1_Index
        HZZ2l2q_resolvedJet2_Index = self.worker.resolvedJet2_Index
        HZZ2l2nu_VBFIndexJet1 = self.worker.HZZ2l2nu_VBFIndexJet1
        HZZ2l2nu_VBFIndexJet2 = self.worker.HZZ2l2nu_VBFIndexJet2
        HZZ2l2nu_minDPhi_METAK4 = self.worker.minDeltaPhi

        # For 2l2nu channel
        HZZ2l2nu_ifVBF = self.worker.HZZ2l2nu_ifVBF
        HZZ2l2qNu_isELE = self.worker.HZZ2l2qNu_isELE
        HZZ2l2nu_isEMuCR = self.worker.HZZ2l2nu_isEMuCR
        HZZ2l2nu_isGammaCR = self.worker.HZZ2l2nu_isGammaCR
        HZZ2l2qNu_cutOppositeChargeFlag = self.worker.HZZ2l2qNu_cutOppositeChargeFlag
        HZZ2l2qNu_nJets = self.worker.HZZ2l2qNu_nJets
        HZZ2l2nu_nExtraLep = self.worker.HZZ2l2nu_nExtraLep
        HZZ2l2nu_nExtraSoftMu = self.worker.HZZ2l2nu_nExtraSoftMu
        HZZ2l2nu_nExtraSoftMuPV = self.worker.HZZ2l2nu_nExtraSoftMuPV
        HZZ2l2nu_passDeltaRl1l2 = self.worker.HZZ2l2nu_passDeltaRl1l2
        HZZ2l2qNu_nTightBtagJets = self.worker.HZZ2l2qNu_nTightBtagJets
        HZZ2l2qNu_nMediumBtagJets = self.worker.HZZ2l2qNu_nMediumBtagJets
        HZZ2l2qNu_nLooseBtagJets = self.worker.HZZ2l2qNu_nLooseBtagJets

        if (foundZZCandidate_4l or foundZZCandidate_2l2q or foundZZCandidate_2l2nu or foundZZCandidate_2l2nu_emuCR or foundGammaCR):
            if self.DEBUG: print("Found ZZ candidate (4l, 2l2q, 2l2nu): ({}, {}, {})".format(foundZZCandidate_4l, foundZZCandidate_2l2q, foundZZCandidate_2l2nu))
            pTL1 = self.worker.pTL1
            etaL1 = self.worker.etaL1
            phiL1 = self.worker.phiL1
            massL1 = self.worker.massL1
            pTL2 = self.worker.pTL2
            etaL2 = self.worker.etaL2
            phiL2 = self.worker.phiL2
            massL2 = self.worker.massL2
            DeltaRl1l2 = self.worker.DeltaRl1l2

            if pTL2>pTL1:
                pTL1, pTL2 = pTL2, pTL1
                etaL1, etaL2 = etaL2, etaL1
                phiL1, phiL2 = phiL2, phiL1
                massL1,massL2 = massL2, massL1

            # Kinematics of Z1: Obtained from pair of leptons with mass closest to Z mass
            pTZ1 = self.worker.Z1.Pt()
            etaZ1 = self.worker.Z1.Eta()
            phiZ1 = self.worker.Z1.Phi()
            massZ1 = self.worker.Z1.M()

            # Kinematics of Z2: Only for 4l and 2l2q channels
            # For 2l2nu channel, Z2 kinematics are obtained from MET
            # For 2l2q channel, Z2 represents the kinamatics of the boosted Z topology
            pTZ2 = self.worker.Z2.Pt()
            etaZ2 = self.worker.Z2.Eta()
            phiZ2 = self.worker.Z2.Phi()
            massZ2 = self.worker.Z2.M()

            pTj1 = self.worker.pTj1
            etaj1 = self.worker.etaj1
            phij1 = self.worker.phij1
            mj1 = self.worker.mj1
            pTj2 = self.worker.pTj2
            etaj2 = self.worker.etaj2
            phij2 = self.worker.phij2
            mj2 = self.worker.mj2

        if (foundZZCandidate_2l2q):
            keepIt = True
            passZZ2l2qSelection = True
            self.passZZ2l2qEvts += 1
            self.CutFlowTable.Fill(4)

            massZ2_2j = self.worker.Z2_2j.M()
            phiZ2_2j = self.worker.Z2_2j.Phi()
            etaZ2_2j = self.worker.Z2_2j.Eta()
            pTZ2_2j = self.worker.Z2_2j.Pt()
            EneZ2_2j = self.worker.Z2_2j.E()

        if (foundZZCandidate_2l2nu):
            keepIt = True
            passZZ2l2nuSelection = True
            self.passZZ2l2nuEvts += 1
            self.CutFlowTable.Fill(5)

            #     FatJet_PNZvsQCD = self.worker.FatJet_PNZvsQCD
            #     self.out.fillBranch("FatJet_PNZvsQCD",FatJet_PNZvsQCD)

        if (foundZZCandidate_2l2nu_emuCR):
            keepIt = True
            passZZ2l2nu_emuCR_Selection = True
            self.passZZ2l2nu_emuCR_Evts += 1
            self.CutFlowTable.Fill(6)

        if (foundGammaCR):
            keepIt = True
            passGammaCR_Selection = True
            self.passGammaCR_Evts += 1
            gamma_pt, gamma_eta, gamma_phi = self.worker.gamma_pt, self.worker.gamma_eta, self.worker.gamma_phi
            for iw, (ptMin, path, weight) in enumerate(self.gammaWindows):
                if gamma_pt >= ptMin:
                    gamma_trigWindowPath = iw
            if gamma_trigWindowPath >= 0:
                _, path, weight = self.gammaWindows[gamma_trigWindowPath]
                gamma_trigWindowFired = fired.get(path, False)
                gamma_prescaleWeight = 1.0 if self.isMC else weight

        if (foundZZCandidate_2l2nu or foundZZCandidate_2l2nu_emuCR or foundGammaCR):
            phiZ2_met = self.worker.Z2_met.Phi()
            pTZ2_met = self.worker.Z2_met.Pt()
            EneZ2_met = self.worker.Z2_met.E()

            # |dphi(Z, MET)|: computed and stored as an output branch; NO cut applied
            # (AN-2016/325 4.4.7 uses |dphi(Z, MET)| > 0.5 - left to a downstream selection).
            HZZ2l2nu_dPhi_ZMET = abs(ROOT.TVector2.Phi_mpi_pi(self.worker.Z1.Phi() - corr_phi))

            HZZ2l2nu_ZZmT = self.worker.MT_2l2nu  # AN-2016/325 Eq. 6
            HZZ2l2nu_ZZpT = self.worker.ZZ_metsystem.Pt()
            

            #Pz_neutrino = self.worker.Pz_neutrino

            # Define TLorentzVector for VBF jets and get dijet mass
            if HZZ2l2nu_VBFIndexJet1>=0 and HZZ2l2nu_VBFIndexJet2>=0:
                VBF_jet1 = ROOT.TLorentzVector()
                VBF_jet2 = ROOT.TLorentzVector()
                VBF_jet1.SetPtEtaPhiM(jets[HZZ2l2nu_VBFIndexJet1].pt, jets[HZZ2l2nu_VBFIndexJet1].eta, jets[HZZ2l2nu_VBFIndexJet1].phi, jets[HZZ2l2nu_VBFIndexJet1].mass)
                VBF_jet2.SetPtEtaPhiM(jets[HZZ2l2nu_VBFIndexJet2].pt, jets[HZZ2l2nu_VBFIndexJet2].eta, jets[HZZ2l2nu_VBFIndexJet2].phi, jets[HZZ2l2nu_VBFIndexJet2].mass)
                VBF_dijet = VBF_jet1 + VBF_jet2
                if self.DEBUG: print("in .py file: VBF_dijet_mass: ", VBF_dijet.M())

                HZZ2l2nu_VBFjet1_pT = jets[HZZ2l2nu_VBFIndexJet1].pt
                HZZ2l2nu_VBFjet1_eta = jets[HZZ2l2nu_VBFIndexJet1].eta
                HZZ2l2nu_VBFjet1_phi = jets[HZZ2l2nu_VBFIndexJet1].phi
                HZZ2l2nu_VBFjet1_mass = jets[HZZ2l2nu_VBFIndexJet1].mass

                HZZ2l2nu_VBFjet2_pT = jets[HZZ2l2nu_VBFIndexJet2].pt
                HZZ2l2nu_VBFjet2_eta = jets[HZZ2l2nu_VBFIndexJet2].eta
                HZZ2l2nu_VBFjet2_phi = jets[HZZ2l2nu_VBFIndexJet2].phi
                HZZ2l2nu_VBFjet2_mass = jets[HZZ2l2nu_VBFIndexJet2].mass

                HZZ2l2nu_VBFdijet_mass = VBF_dijet.M()
                HZZ2l2nu_VBFdijet_pT = VBF_dijet.Pt()
                HZZ2l2nu_VBFdijet_E = VBF_dijet.E()

                HZZ2l2nu_VBFdEta_jj = abs(VBF_jet1.Eta() - VBF_jet2.Eta())
                HZZ2l2nu_VBFdPhi_jj = abs(VBF_jet1.DeltaPhi(VBF_jet2))
                HZZ2l2nu_VBFdR_jj = VBF_jet1.DeltaR(VBF_jet2)

            # print("inside 2l2nu loop")
        # self.out.fillBranch("phiZ2_met",phiZ2_met)
        # self.out.fillBranch("pTZ2_met",pTZ2_met)
        # self.out.fillBranch("EneZ2_met",EneZ2_met)

        #if (foundZZCandidate_2l2nu_emuCR):
            #keepIt = True
            #passZZ2l2nu_emuCR_Selection = True
            #self.passZZ2l2nu_emuCR_Evts += 1


        if (foundZZCandidate_4l):
            keepIt = True
            self.passZZ4lEvts += 1
            self.CutFlowTable.Fill(3)
            passZZ4lSelection = True
            # print("Inside 4l loop: ",passZZ2l2qSelection)
            D_CP = self.worker.D_CP
            D_0m = self.worker.D_0m
            D_0hp = self.worker.D_0hp
            D_int = self.worker.D_int
            D_L1 = self.worker.D_L1
            D_L1Zg = self.worker.D_L1Zg

            pTL3 = self.worker.pTL3
            etaL3 = self.worker.etaL3
            phiL3 = self.worker.phiL3
            massL3 = self.worker.massL3
            pTL4 = self.worker.pTL4
            etaL4 = self.worker.etaL4
            phiL4 = self.worker.phiL4
            massL4 = self.worker.massL4
            

            if pTL4>pTL3:
                pTL3, pTL4 = pTL4, pTL3
                etaL3, etaL4 = etaL4, etaL3
                phiL3, phiL4 = phiL4, phiL3
                massL3, massL4 = massL4, massL3

            pT4l = self.worker.ZZsystem.Pt()
            eta4l = self.worker.ZZsystem.Eta()
            phi4l = self.worker.ZZsystem.Phi()
            mass4l = self.worker.ZZsystem.M()
            if self.worker.isFSR==False:
                pT4l = self.worker.ZZsystemnofsr.Pt()
                eta4l = self.worker.ZZsystemnofsr.Eta()
                phi4l = self.worker.ZZsystemnofsr.Phi()
                mass4l = self.worker.ZZsystemnofsr.M()

        if self.DEBUG:
            print("(found candidates: 2l2q, 2l2nu, 4l): ({:1}, {:1}, {:1}), pTL1: {:7.3f}, pTL2: {:7.3f}, pTZ1: {:7.3f}, pTZ2: {:7.3f}, pTZ2_2j: {:7.3f}, pTZ2_met: {:7.3f}".format(foundZZCandidate_2l2q, foundZZCandidate_2l2nu, foundZZCandidate_4l, pTL1, pTL2, pTZ1, pTZ2, pTZ2_2j, pTZ2_met))
            if (foundZZCandidate_2l2q):
                print("==> pTL1: {}, \t pTL2: {}".format(pTL1, pTL2))

        # Fill the branches with the Trigger information for each channel
        for TriggerChannel in self.cfg['TriggerChannels']:
            self.out.fillBranch(TriggerChannel, TriggerMap[TriggerChannel])
        self.out.fillBranch("phiZ2_met",phiZ2_met)
        self.out.fillBranch("pTZ2_met",pTZ2_met)
        self.out.fillBranch("EneZ2_met",EneZ2_met)
        self.out.fillBranch("HZZ2l2nu_ZZmT", HZZ2l2nu_ZZmT)
        self.out.fillBranch("HZZ2l2nu_ZZpT", HZZ2l2nu_ZZpT)
        self.out.fillBranch("HZZ2l2nu_minDPhi_METAK4", HZZ2l2nu_minDPhi_METAK4)
        self.out.fillBranch("HZZ2l2nu_dPhi_ZMET", HZZ2l2nu_dPhi_ZMET)
        #self.out.fillBranch("Pz_neutrino", Pz_neutrino)

        self.out.fillBranch("HZZ2l2nu_VBFIndexJet1", HZZ2l2nu_VBFIndexJet1)
        self.out.fillBranch("HZZ2l2nu_VBFIndexJet2", HZZ2l2nu_VBFIndexJet2)

        self.out.fillBranch("HZZ2l2nu_VBFjet1_pT", HZZ2l2nu_VBFjet1_pT)
        self.out.fillBranch("HZZ2l2nu_VBFjet1_eta", HZZ2l2nu_VBFjet1_eta)
        self.out.fillBranch("HZZ2l2nu_VBFjet1_phi", HZZ2l2nu_VBFjet1_phi)
        self.out.fillBranch("HZZ2l2nu_VBFjet1_mass", HZZ2l2nu_VBFjet1_mass)

        self.out.fillBranch("HZZ2l2nu_VBFjet2_pT", HZZ2l2nu_VBFjet2_pT)
        self.out.fillBranch("HZZ2l2nu_VBFjet2_eta", HZZ2l2nu_VBFjet2_eta)
        self.out.fillBranch("HZZ2l2nu_VBFjet2_phi", HZZ2l2nu_VBFjet2_phi)
        self.out.fillBranch("HZZ2l2nu_VBFjet2_mass", HZZ2l2nu_VBFjet2_mass)

        self.out.fillBranch("HZZ2l2nu_VBFdijet_mass", HZZ2l2nu_VBFdijet_mass)
        self.out.fillBranch("HZZ2l2nu_VBFdijet_pT", HZZ2l2nu_VBFdijet_pT)
        self.out.fillBranch("HZZ2l2nu_VBFdijet_E", HZZ2l2nu_VBFdijet_E)
        self.out.fillBranch("HZZ2l2nu_VBFdEta_jj", HZZ2l2nu_VBFdEta_jj)
        self.out.fillBranch("HZZ2l2nu_VBFdPhi_jj", HZZ2l2nu_VBFdPhi_jj)
        self.out.fillBranch("HZZ2l2nu_VBFdR_jj", HZZ2l2nu_VBFdR_jj)

        self.out.fillBranch("massZ2_2j",massZ2_2j)
        self.out.fillBranch("phiZ2_2j",phiZ2_2j)
        self.out.fillBranch("etaZ2_2j",etaZ2_2j)
        self.out.fillBranch("pTZ2_2j",pTZ2_2j)
        self.out.fillBranch("EneZ2_2j",EneZ2_2j)

        self.out.fillBranch("massL1",massL1)
        self.out.fillBranch("pTL1",pTL1)
        self.out.fillBranch("etaL1",etaL1)
        self.out.fillBranch("phiL1",phiL1)
        self.out.fillBranch("massL2",massL2)
        self.out.fillBranch("pTL2",pTL2)
        self.out.fillBranch("etaL2",etaL2)
        self.out.fillBranch("phiL2",phiL2)
        self.out.fillBranch("DeltaRl1l2",DeltaRl1l2)
        self.out.fillBranch("HZZ2l2nu_passDeltaRl1l2",HZZ2l2nu_passDeltaRl1l2)

        self.out.fillBranch("pTZ1",pTZ1)
        self.out.fillBranch("etaZ1",etaZ1)
        self.out.fillBranch("phiZ1",phiZ1)
        self.out.fillBranch("massZ1",massZ1)
        self.out.fillBranch("pTZ2",pTZ2)
        self.out.fillBranch("etaZ2",etaZ2)
        self.out.fillBranch("phiZ2",phiZ2)
        self.out.fillBranch("massZ2",massZ2)
        self.out.fillBranch("pT_MET",corr_pt) #for v15
        self.out.fillBranch("phi_MET",corr_phi) #for v15

        self.out.fillBranch("overallEventWeight",             overallEventWeight)
        self.out.fillBranch("overallEventWeight_puUp",        overallEventWeight_puUp)
        self.out.fillBranch("overallEventWeight_puDown",      overallEventWeight_puDown)
        self.out.fillBranch("overallEventWeight_prefireUp",   overallEventWeight_prefireUp)
        self.out.fillBranch("overallEventWeight_prefireDown", overallEventWeight_prefireDown)

        self.out.fillBranch("mass4l",mass4l)
        self.out.fillBranch("pT4l",pT4l)
        self.out.fillBranch("eta4l",eta4l)
        self.out.fillBranch("phi4l",phi4l)
        self.out.fillBranch("D_CP",D_CP)
        self.out.fillBranch("D_0m",D_0m)
        self.out.fillBranch("D_0hp",D_0hp)
        self.out.fillBranch("D_int",D_int)
        self.out.fillBranch("D_L1",D_L1)
        self.out.fillBranch("D_L1Zg",D_L1Zg)

        self.out.fillBranch("massL3",massL3)
        self.out.fillBranch("pTL3",pTL3)
        self.out.fillBranch("etaL3",etaL3)
        self.out.fillBranch("phiL3",phiL3)
        self.out.fillBranch("massL4",massL4)
        self.out.fillBranch("pTL4",pTL4)
        self.out.fillBranch("etaL4",etaL4)
        self.out.fillBranch("phiL4",phiL4)

        self.out.fillBranch("mj1",mj1)
        self.out.fillBranch("pTj1",pTj1)
        self.out.fillBranch("etaj1",etaj1)
        self.out.fillBranch("phij1",phij1)
        self.out.fillBranch("mj2",mj2)
        self.out.fillBranch("pTj2",pTj2)
        self.out.fillBranch("etaj2",etaj2)
        self.out.fillBranch("phij2",phij2)

        self.out.fillBranch("passZZ4lSelection",passZZ4lSelection)
        self.out.fillBranch("passZZ2l2qSelection",passZZ2l2qSelection)
        self.out.fillBranch("passZZ2l2nuSelection",passZZ2l2nuSelection)
        self.out.fillBranch("passZZ2l2nu_emuCR_Selection",passZZ2l2nu_emuCR_Selection)
        self.out.fillBranch("isBoosted2l2q",isBoosted2l2q)

        self.out.fillBranch("HZZ2l2q_boostedJet_PNScore", HZZ2l2q_boostedJet_PNScore)
        self.out.fillBranch("HZZ2l2q_boostedJet_Index", HZZ2l2q_boostedJet_Index)

        self.out.fillBranch("HZZ2l2q_resolvedJet1_Index",HZZ2l2q_resolvedJet1_Index)
        self.out.fillBranch("HZZ2l2q_resolvedJet2_Index",HZZ2l2q_resolvedJet2_Index)

        self.out.fillBranch("HZZ2l2nu_ifVBF",HZZ2l2nu_ifVBF)
        self.out.fillBranch("HZZ2l2qNu_isELE",HZZ2l2qNu_isELE)
        self.out.fillBranch("HZZ2l2nu_isEMuCR",HZZ2l2nu_isEMuCR)
        self.out.fillBranch("passGammaCR_Selection",passGammaCR_Selection)
        self.out.fillBranch("HZZ2l2nu_isGammaCR",HZZ2l2nu_isGammaCR)
        self.out.fillBranch("gamma_pt",gamma_pt)
        self.out.fillBranch("gamma_eta",gamma_eta)
        self.out.fillBranch("gamma_phi",gamma_phi)
        self.out.fillBranch("gamma_trigWindowFired",gamma_trigWindowFired)
        self.out.fillBranch("gamma_trigWindowPath",gamma_trigWindowPath)
        self.out.fillBranch("gamma_prescaleWeight",gamma_prescaleWeight)
        self.out.fillBranch("HZZ2l2qNu_cutOppositeChargeFlag",HZZ2l2qNu_cutOppositeChargeFlag)
        self.out.fillBranch("HZZ2l2qNu_nJets",HZZ2l2qNu_nJets)
        self.out.fillBranch("HZZ2l2nu_nExtraLep",HZZ2l2nu_nExtraLep)
        self.out.fillBranch("HZZ2l2nu_nExtraSoftMu",HZZ2l2nu_nExtraSoftMu)
        self.out.fillBranch("HZZ2l2nu_nExtraSoftMuPV",HZZ2l2nu_nExtraSoftMuPV)
        self.out.fillBranch("HZZ2l2qNu_nTightBtagJets",HZZ2l2qNu_nTightBtagJets)
        self.out.fillBranch("HZZ2l2qNu_nMediumBtagJets",HZZ2l2qNu_nMediumBtagJets)
        self.out.fillBranch("HZZ2l2qNu_nLooseBtagJets",HZZ2l2qNu_nLooseBtagJets)

        # NOTE: overallEventWeight (= genWeight * puWeight * L1PreFiringWeight_Nom)
        # and its PU / prefiring up-down variations are filled above.
        # STILL MISSING (need a POG payload + verified WP key, or a per-variation
        # re-run of the C++ selection): lepton reco/ID/iso/trigger SFs, b-tag SF,
        # electron energy scale/smearing, JES/JER/unclustered-MET shape variations,
        # Rochester muon-scale propagation.

        return keepIt

# define modules using the syntax 'name = lambda : constructor' to avoid
# having them loaded when not needed

# H4LCppModulie() = lambda: HZZAnalysisCppProducer(year)
