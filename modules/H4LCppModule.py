from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
import ROOT
import yaml
import os
from collections import OrderedDict

from modules.Helper import *
from modules.METFilters import passFilters

ROOT.PyConfig.IgnoreCommandLineOptions = True


class HZZAnalysisCppProducer(Module):

    def __init__(self, year, cfgFile, isMC, isFSR, channels, DEBUG=False):
        self.loadLibraries()
        self.year = year
        self.isMC = isMC
        self.channels = channels # choices=["all", "4l", "2l2q", "2l2v"],
        self.DEBUG = DEBUG
        self.cfgFile = cfgFile
        self.cfg = self._load_config(cfgFile)
        self.worker = ROOT.H4LTools(int(self.year), self.isMC, self.DEBUG)
        self.genworker = ROOT.GenAnalysis(self.DEBUG)
        self._initialize_worker(self.cfg)
        self.worker.isFSR = isFSR
        self._initialize_counters()
        if str(self.year) == '2022':
            self.PUweight_list = self.GetPUWeight()

        # Alternatively, for dynamic worker attributes
        self.dynamicCuts_4l = ["cut4e", "cutghost4e", "cutLepPt4e", "cutQCD4e", "cutZZ4e", "cutm4l4e",
                          "cut4mu", "cutghost4mu", "cutLepPt4mu", "cutQCD4mu", "cutZZ4mu", "cutm4l4mu",
                          "cut2e2mu", "cutghost2e2mu", "cutLepPt2e2mu", "cutQCD2e2mu", "cutZZ2e2mu",
                          "cutm4l2e2mu"]
        self.dynamicCuts_2l2q = ["HZZ2l2qNu_cut2l", "HZZ2l2qNu_cutOppositeCharge", "HZZ2l2qNu_cutpTl1l2",
                             "HZZ2l2qNu_cutETAl1l2", "HZZ2l2qNu_cutmZ1Window", "HZZ2l2qNu_cutZ1Pt",
                             "cut2l1J", "cut2l2j", "cut2l1Jor2j"]
        self.dynamicCuts_2l2nu = ["HZZ2l2qNu_cut2l", "HZZ2l2qNu_cutOppositeCharge", "HZZ2l2qNu_cutpTl1l2",
                             "HZZ2l2qNu_cutETAl1l2", "HZZ2l2qNu_cutmZ1Window", "HZZ2l2qNu_cutZ1Pt",
                             "HZZ2l2nu_cutbtag", "HZZ2l2nu_cutdPhiJetMET", "HZZ2l2nu_cutMETgT100"]
        self.dynamicCuts_2l2nu_emu_CR = ["HZZemuCR_cut2l", "HZZemuCR_cutpTl1l2",
                             "HZZemuCR_cutETAl1l2", "HZZemuCR_cutmZ1Window",
                             "HZZemuCR_cutZ1Pt", "HZZ2l2nu_cutdPhiJetMET", "HZZ2l2nu_cutMETgT100"]


        hist_nBins = 7 + len(self.dynamicCuts_4l) + len(self.dynamicCuts_2l2q) + len(self.dynamicCuts_2l2nu) + len(self.dynamicCuts_2l2nu_emu_CR)
        print("hist_nBins: ", hist_nBins)
        self.CutFlowTable =  ROOT.TH1F('cutFlow','cutFlow',hist_nBins, 0, hist_nBins)
        self.CutFlowTable.GetXaxis().SetBinLabel(1, "Total")
        self.CutFlowTable.GetXaxis().SetBinLabel(2, "PassTrig")
        self.CutFlowTable.GetXaxis().SetBinLabel(3, "PassMETFilters")
        self.CutFlowTable.GetXaxis().SetBinLabel(4, "PassZZSelection")
        self.CutFlowTable.GetXaxis().SetBinLabel(5, "PassZZ2l2qSelection")
        self.CutFlowTable.GetXaxis().SetBinLabel(6, "PassZZ2l2nuSelection")
        self.CutFlowTable.GetXaxis().SetBinLabel(7, "PassZZ2l2nu_emuCR_Selection")

        for idx, cut in enumerate(self.dynamicCuts_4l):
            self.CutFlowTable.GetXaxis().SetBinLabel(8 + idx, cut)
        for idx, cut in enumerate(self.dynamicCuts_2l2q):
            self.CutFlowTable.GetXaxis().SetBinLabel(8 + len(self.dynamicCuts_4l) + idx, cut)
        for idx, cut in enumerate(self.dynamicCuts_2l2nu):
            self.CutFlowTable.GetXaxis().SetBinLabel(8 + len(self.dynamicCuts_4l) + len(self.dynamicCuts_2l2q) + idx, cut)
        for idx, cut in enumerate(self.dynamicCuts_2l2nu_emu_CR):
            self.CutFlowTable.GetXaxis().SetBinLabel(8 + len(self.dynamicCuts_4l) + len(self.dynamicCuts_2l2q) + len(self.dynamicCuts_2l2nu) + idx, cut)

    def GetPUWeight(self):
        PUinput_file = ROOT.TFile.Open(self.PUweightfile)
        PUinput_hist = PUinput_file.Get(self.PUweighthisto)
        PUweight_list = []
        for i in range(1, PUinput_hist.GetNbinsX() + 1):
            PUweight_list.append(PUinput_hist.GetBinContent(i))
        PUinput_file.Close()
        return PUweight_list

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

        # Load the C++ GEN module
        if "/GenAnalysis_cc.so" not in ROOT.gSystem.GetLibraries():
            print("Load GenAnalysis C++ module")
            if base_path:
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/GenAnalysis.cc+O" % base_path)
            else:
                base_path = "$CMSSW_BASE//src/PhysicsTools/NanoAODTools"
                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gROOT.ProcessLine(
                    ".L %s/include/GenAnalysis.h" % base_path)

    def _load_config(self, cfgFile):
        with open(cfgFile, 'r') as ymlfile:
            return yaml.safe_load(ymlfile)

    def _initialize_worker(self, cfg):
        print("Initializing worker")
        if self.DEBUG:
            print("Config file loaded: ", cfg)
        if str(self.year) == '2022':
            self.PUweightfile = cfg["outputdataNPV"]
            self.PUweighthisto = cfg["PUweightHistoName"]
        self.worker.InitializeElecut(*self._get_nested_values(cfg['Electron'], [
            'pTcut', 'Etacut', 'Sip3dcut', 'Loosedxycut', 'Loosedzcut',
            'Isocut', ['BDTWP', 'LowEta', 'LowPT'], ['BDTWP', 'MedEta', 'LowPT'],
            ['BDTWP', 'HighEta', 'LowPT'], ['BDTWP', 'LowEta', 'HighPT'],
            ['BDTWP', 'MedEta', 'HighPT'], ['BDTWP', 'HighEta', 'HighPT']
        ]))

        self.worker.InitializeMucut(*self._get_nested_values(cfg['Muon'], [
            'pTcut', 'Etacut', 'Sip3dcut', 'Loosedxycut', 'Loosedzcut', 'Isocut',
            'Tightdxycut', 'Tightdzcut', 'TightTrackerLayercut', 'TightpTErrorcut',
            'HighPtBound'
        ]))

        self.worker.InitializeFsrPhotonCut(*self._get_nested_values(cfg['FsrPhoton'], [
            'pTcut', 'Etacut', 'Isocut', 'dRlcut', 'dRlOverPtcut'
        ]))

        self.worker.InitializeJetcut(*self._get_nested_values(cfg['Jet'], ['pTcut', 'Etacut']))

        self.worker.InitializeEvtCut(*self._get_nested_values(cfg, ['MZ1cut', 'MZZcut',
                                                                    ['Higgscut', 'down'], ['Higgscut', 'up'],
                                                                    'Zmass', ['MZcut', 'down'], ['MZcut', 'up'],
                                                                    ['Jet','deepJet_btag','Loose'], ['Jet','deepJet_btag','Medium'], ['Jet','deepJet_btag','Tight']]))

        if self.channels == "all" or self.channels == "2l2q":
            self.worker.InitializeHZZ2l2qCut(*self._get_nested_values(cfg['HZZ2l2q'],
                                                                  ['Leading_Lep_pT', 'SubLeading_Lep_pT', 'Lep_eta',
                                                                   ['MZLepcut', 'down'], ['MZLepcut', 'up']]))

        if self.channels == "all" or self.channels == "2l2v":
            self.worker.InitializeHZZ2l2nuCut(*self._get_nested_values(cfg['HZZ2l2nu'],
                                                                   ['Leading_Lep_pT', 'SubLeading_Lep_pT', 'Lep_eta', 'Pt_ll',
                                                                    'M_ll_Window', 'dPhi_jetMET', ['MZLepcut', 'down'], ['MZLepcut', 'up']]))

        print("Worker initialized")


    def _get_nested_values(self, dictionary, keys):
        values = []
        for key in keys:
            if isinstance(key, list):
                sub_dict = dictionary
                for sub_key in key:
                    sub_dict = sub_dict.get(sub_key, {})
                values.append(sub_dict if sub_dict else 'N/A')
            else:
                values.append(dictionary.get(key, 'N/A'))
        return values

    def _initialize_counters(self):
        self.passAllEvts = 0
        self.passtrigEvts = 0
        self.passMETFilters = 0
        self.passZZ4lEvts = 0
        self.passZZ2l2qEvts = 0
        self.passZZ2l2nuEvts = 0
        self.passZZ2l2nu_emuCR_Evts = 0

    def beginJob(self):
        pass

    def endJob(self):
        print("PassTrig: "+str(self.passtrigEvts)+" Events")
        print("Pass4eCut: "+str(self.worker.cut4e)+" Events")
        print("Pass4eGhostRemoval: "+str(self.worker.cutghost4e)+" Events")
        print("Pass4eLepPtCut: "+str(self.worker.cutLepPt4e)+" Events")
        print("Pass4eQCDSupress: "+str(self.worker.cutQCD4e)+" Events")
        print("PassmZ1mZ2Cut_4e: "+str(self.worker.cutZZ4e)+" Events")
        print("Passm4l_105_160_Cut_4e: "+str(self.worker.cutm4l4e)+" Events")
        print("Pass4muCut: "+str(self.worker.cut4mu)+" Events")
        print("Pass4muGhostRemoval: "+str(self.worker.cutghost4mu)+" Events")
        print("Pass4muLepPtCut: "+str(self.worker.cutLepPt4mu)+" Events")
        print("Pass4muQCDSupress: "+str(self.worker.cutQCD4mu)+" Events")
        print("PassmZ1mZ2Cut_4mu: "+str(self.worker.cutZZ4mu)+" Events")
        print("Passm4l_105_160_Cut_4mu: "+str(self.worker.cutm4l4mu)+" Events")
        print("Pass2e2muCut: "+str(self.worker.cut2e2mu)+" Events")
        print("Pass2e2muGhostRemoval: "+str(self.worker.cutghost2e2mu)+" Events")
        print("Pass2e2muLepPtCut: "+str(self.worker.cutLepPt2e2mu)+" Events")
        print("Pass2e2muQCDSupress: "+str(self.worker.cutQCD2e2mu)+" Events")
        print("PassmZ1mZ2Cut_2e2mu: "+str(self.worker.cutZZ2e2mu)+" Events")
        print("Passm4l_105_160_Cut_2e2mu: "+str(self.worker.cutm4l2e2mu)+" Events")
        print("PassZZSelection: "+str(self.passZZEvts)+" Events")
        if self.isMC:
            print("PassGEN4eCut: "+str(self.genworker.nGEN4e)+" Events")
            print("PassGEN4eZ1Cut: "+str(self.genworker.nGEN4epassZ1)+" Events")
            print("PassGEN4efidCut: "+str(self.genworker.nGEN4epassFid)+" Events")
            print("PassGEN2e2muCut: "+str(self.genworker.nGEN2e2mu)+" Events")
            print("PassGEN2e2muZ1Cut: "+str(self.genworker.nGEN2e2mupassZ1)+" Events")
            print("PassGEN2e2mufidCut: "+str(self.genworker.nGEN2e2mupassFid)+" Events")
            print("PassGEN4muCut: "+str(self.genworker.nGEN4mu)+" Events")
            print("PassGEN4muZ1Cut: "+str(self.genworker.nGEN4mupassZ1)+" Events")
            print("PassGEN4mufidCut: "+str(self.genworker.nGEN4mupassFid)+" Events")

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree)  # initReaders must be called in beginFile
        self.out = wrappedOutputTree

        # Boolean branches for Trigger channels
        for TriggerChannel in self.cfg['TriggerChannels']:
            self.out.branch(TriggerChannel, "O")
        # boolean branches for 4l, 2l2q, 2l2nu channels
        self.out.branch("foundZZCandidate_4l",  "O")
        self.out.branch("foundZZCandidate_2l2q",  "O")
        self.out.branch("foundZZCandidate_2l2nu",  "O")
        self.out.branch("passZZ2l2nu_emuCR_Selection", "O")
        self.out.branch("isBoosted2l2q",  "O")
        self.out.branch("HZZ2l2nu_ifVBF", "O")
        self.out.branch("HZZ2l2qNu_isELE", "O")
        self.out.branch("HZZ2l2qNu_cutOppositeChargeFlag", "O")
        self.out.branch("HZZ2l2nu_isEMuCR", "O")

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

        # Branches for 4l channel: ZZ kinematics
        self.out.branch("mass4l",  "F")
        self.out.branch("GENmass4l",  "F")
        self.out.branch("mass4e",  "F")
        self.out.branch("mass4mu",  "F")
        self.out.branch("mass2e2mu",  "F")
        self.out.branch("pT4l",  "F")
        self.out.branch("GENpT4l",  "F")
        self.out.branch("rapidity4l",  "F")
        self.out.branch("njets_pt30_eta4p7", "I")
        self.out.branch("finalState", "I")
        self.out.branch("GENnjets_pt30_eta4p7", "I")
        self.out.branch("GENrapidity4l",  "F")
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


        self.out.branch("mj2",  "F")
        self.out.branch("EvtNum",  "I")
        self.out.branch("Weight",  "F")
        self.out.branch("pileupWeight",  "F")
        self.out.branch("dataMCWeight_new",  "F")
        self.out.branch("prefiringWeight",  "F")
        self.out.branch("passedTrig",  "O")
        self.out.branch("passedFullSelection",  "O")
        self.out.branch("passedZ4lSelection",  "O")
        self.out.branch("passedZ4lZ1LSelection",  "O")
        self.out.branch("passedZ4lZXCRSelection",  "O")
        self.out.branch("passedZXCRSelection",  "O")
        self.out.branch("passedFiducialSelection",  "O")
        GENHlepNum = 4
        GENZNum = 2
        self.out.branch("GENlep_MomId",  "I", lenVar = "nGenPart")
        self.out.branch("GENlep_MomMomId",  "I", lenVar = "nGenPart")
        self.out.branch("GENZ_MomId",  "I", lenVar = "nVECZ")
        self.out.branch("GENZ_DaughtersId",  "I", lenVar = "GENZNum")
        self.out.branch("GENlep_Hindex",  "I", lenVar = "GENHlepNum")
        self.out.branch("lep_Hindex",  "I", lenVar = "GENHlepNum")
        self.out.branch("GENlep_id",  "I", lenVar = "nGENLeptons")
        self.out.branch("lep_genindex",  "I", lenVar = "Lepointer")
        # FSR branches for leptons
        self.out.branch("Electron_Fsr_pt",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_eta",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_phi",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Muon_Fsr_pt",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_eta",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_phi",  "F", lenVar = "nMuon_Fsr")

        # Branches for 2l2q channel
        self.out.branch("massZ2_2j",  "F")
        self.out.branch("phiZ2_2j",  "F")
        self.out.branch("etaZ2_2j",  "F")
        self.out.branch("pTZ2_2j",  "F")
        self.out.branch("EneZ2_2j",  "F")
        self.out.branch("phiZ2_met",  "F")
        self.out.branch("pTZ2_met",  "F")
        self.out.branch("EneZ2_met",  "F")
        self.out.branch("MT_2l2nu",  "F")

        # Branches for 2l2nu channel: ZZ kinematics
        self.out.branch("HZZ2l2nu_ZZmT",  "F")
        self.out.branch("HZZ2l2nu_ZZpT",  "F")
        #self.out.branch("Pz_neutrino", "F")

        # Branches for 2l2nu channel: VBF jets and dijet kinematics
        self.out.branch("HZZ2l2qNu_nJets", "I")

        self.out.branch("HZZ2l2qNu_nTightBtagJets", "I")
        self.out.branch("HZZ2l2qNu_nMediumBtagJets", "I")
        self.out.branch("HZZ2l2qNu_nLooseBtagJets", "I")

        self.out.branch("HZZ2l2nu_VBFIndexJet1",  "I")
        self.out.branch("HZZ2l2nu_VBFIndexJet2",  "I")
        self.out.branch("HZZ2l2nu_minDPhi_METAK4",  "F")

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


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
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
            self.genworker.Initialize()
            self.genworker.SetObjectNumGen(event.nGenPart, event.nGenJet)
        keepIt = False
        self.passAllEvts += 1
        self.CutFlowTable.Fill(0)
        Lepointer = 0
        EvtNum = 0
        Weight = 1

        foundZZCandidate_4l = False    # for 4l
        foundZZCandidate_2l2q = False # for 2l2q
        foundZZCandidate_2l2nu = False # for 2l2nu
        foundZZCandidate_2l2nu_emuCR = False   # for 2l2nu emu control region
        passZZ2l2nu_emuCR_Selection = False    # for 2l2nu emu control region
        HZZ2l2nu_ifVBF = False
        HZZ2l2qNu_isELE = False
        HZZ2l2nu_isEMuCR = False
        HZZ2l2qNu_cutOppositeChargeFlag = False
        isBoosted2l2q = False

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
        prefiringWeight = 1
        dataMCWeight_new = 1
        pileupWeight = 1
        mass4e=0
        mass2e2mu=0
        mass4mu=0
        finalState=-1
        GENmass4l = -99
        GENpT4l = -99
        nVECZ = 2
        GENrapidity4l = -99
        GENnjets_pt30_eta4p7 = -1
        nGENLeptons = 0
        nGenPart = 0
        pTZ1 = -99
        etaZ1 = -99
        phiZ1 = -99
        massZ1 = 0
        pTZ2 = -99
        etaZ2 = -99
        phiZ2 = -99
        massZ2 = -99
        pT4l = -99
        eta4l = -99
        phi4l = -99
        mass4l = 0
        rapidity4l = -99

        massZ2_2j = -999.
        phiZ2_2j = -999.
        etaZ2_2j = -999.
        pTZ2_2j = -999.
        EneZ2_2j = -999.
        phiZ2_met = -999.
        pTZ2_met = -999.
        EneZ2_met = -999.
        MT_2l2nu = -999.
        HZZ2l2nu_minDPhi_METAK4 = 999.0

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
        #Pz_neutrino = -999.

        TriggerMap = {}
        passedTrig = False
        for TriggerChannel in self.cfg['TriggerChannels']:
            TriggerMap[TriggerChannel] = PassTrig(event, self.cfg, TriggerChannel)

        # If any of the trigger channel from TriggerMap passes, then the event is kept else return keepIt
        for value in TriggerMap.values():
            if value:
                passedTrig = True
                break
        if not passedTrig:
            return keepIt
        self.passtrigEvts += 1
        self.CutFlowTable.Fill(1)

        if passFilters(event, int(self.year)):
            self.passMETFilters += 1
            self.CutFlowTable.Fill(2)
        else:
            return keepIt

        if(isMC and self.year == 2022):
            pileupWeight = self.PUweight_list[event.Pileup_nPU]
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        # Photons = Collection(event, "Photon")
        jets = Collection(event, "Jet")
        FatJets = Collection(event, "FatJet")
        met = Object(event, "MET", None)

        # for photon in Photons:
        #     # Keep photons if pT > 55, |eta| < 2.5 and skip the transition region of barrel and endcap
        #     if photon.pt > 55 and abs(photon.eta) < 2.5 and not (1.4442 < abs(photon.eta) < 1.566):
        #         return keepIt

        if isMC:
            nGenPart = event.nGenPart
            genparts = Collection(event, "GenPart")
            genjets = Collection(event, "GenJet")
            for xj in genjets:
                self.genworker.SetGenJets(xj.pt,xj.eta,xj.phi,xj.mass)
            for xg in genparts:
                self.worker.SetGenParts(xg.pt)
                self.genworker.SetGenParts(xg.pt,xg.eta,xg.phi,xg.mass,xg.pdgId,xg.status,xg.statusFlags,xg.genPartIdxMother)
            for xm in muons:
                self.worker.SetMuonsGen(xm.genPartIdx)
            for xe in electrons:
                self.worker.SetElectronsGen(xe.genPartIdx)

        for xe in electrons:
            if str(self.year) == "2022":
                self.worker.SetElectrons(xe.pt, xe.eta, xe.phi, xe.mass, xe.dxy,
                                      xe.dz, xe.sip3d, xe.mvaHZZIso, xe.mvaIso_WP90, xe.pdgId, xe.pfRelIso03_all)
            else:
                self.worker.SetElectrons(xe.pt, xe.eta, xe.phi, xe.mass, xe.dxy,
                                      xe.dz, xe.sip3d, xe.mvaFall17V2Iso, xe.mvaFall17V2Iso_WP90, xe.pdgId, xe.pfRelIso03_all)
            if self.DEBUG:
                print("Electrons: pT, eta: {}, {}".format(xe.pt, xe.eta))

        for xm in muons:
            self.worker.SetMuons(xm.corrected_pt, xm.eta, xm.phi, xm.mass, xm.isGlobal, xm.isTracker, xm.mediumId,
                                xm.dxy, xm.dz, xm.sip3d, xm.ptErr, xm.nTrackerLayers, xm.isPFcand,
                                 xm.pdgId, xm.charge, xm.pfRelIso03_all)
            if self.DEBUG:
                print("Muons: pT, eta: {}, {}".format(xm.corrected_pt, xm.eta))

        for xf in fsrPhotons:
            if self.DEBUG:
                print("FsrPhotons: pT, eta: {}, {}".format(xf.pt, xf.eta))
                print("Year: ", self.year)
            if str(self.year) == '2022':
                self.worker.SetFsrPhotons(xf.dROverEt2, xf.eta, xf.phi, xf.pt, xf.relIso03, xf.electronIdx, xf.muonIdx)
            else:
                self.worker.SetFsrPhotons(xf.dROverEt2, xf.eta, xf.phi, xf.pt, xf.relIso03, 0, 0) # FIXME: need to update the FSR photon index

        for xj in jets:
            if str(self.year) == '2022':
                self.worker.SetJets(xj.pt, xj.eta, xj.phi, xj.mass, xj.jetId, 0.8, 7) # FIXME: need to update the btagDeepFlavB and puId currently these branches are not available
            else:
                self.worker.SetJets(xj.pt, xj.eta, xj.phi, xj.mass, xj.jetId, xj.btagDeepFlavB, xj.puId)

        for xj in FatJets:
            ZvsQCD = xj.particleNetWithMass_ZvsQCD if str(self.year) == '2022' else xj.particleNet_ZvsQCD
            self.worker.SetFatJets(xj.pt, xj.eta, xj.phi, xj.msoftdrop, xj.jetId, xj.btagDeepB, ZvsQCD)

        self.worker.SetMET(met.pt,met.phi,met.sumEt)
        self.worker.BatchFsrRecovery_Run3()

        self.worker.LeptonSelection()
        if self.DEBUG:
            print("Lepton selection done")

        if isMC:
            self.genworker.SetGenVariables()
            GENmass4l = self.genworker.GENmass4l
            GENpT4l = self.genworker.GENpT4l
            GENrapidity4l = self.genworker.GENrapidity4l
            GENnjets_pt30_eta4p7 = self.genworker.GENnjets_pt30_eta4p7
            nGENLeptons = self.genworker.nGENLeptons
            passedFiducialSelection = self.genworker.passedFiducialSelection

        Electron_Fsr_pt_vec = self.worker.ElectronFsrPt()
        Electron_Fsr_eta_vec = self.worker.ElectronFsrEta()
        Electron_Fsr_phi_vec = self.worker.ElectronFsrPhi()
        Muon_Fsr_pt_vec = self.worker.MuonFsrPt()
        Muon_Fsr_eta_vec = self.worker.MuonFsrEta()
        Muon_Fsr_phi_vec = self.worker.MuonFsrPhi()

        Electron_Fsr_pt = []
        Electron_Fsr_eta = []
        Electron_Fsr_phi = []
        Muon_Fsr_pt = []
        Muon_Fsr_eta = []
        Muon_Fsr_phi = []

        if len(Electron_Fsr_pt_vec)>0:
            for i in range(len(Electron_Fsr_pt_vec)):
                Electron_Fsr_pt.append(Electron_Fsr_pt_vec[i])
                Electron_Fsr_eta.append(Electron_Fsr_eta_vec[i])
                Electron_Fsr_phi.append(Electron_Fsr_phi_vec[i])
        if len(Muon_Fsr_pt_vec)>0:
            for i in range(len(Muon_Fsr_pt_vec)):
                Muon_Fsr_pt.append(Muon_Fsr_pt_vec[i])
                Muon_Fsr_eta.append(Muon_Fsr_eta_vec[i])
                Muon_Fsr_phi.append(Muon_Fsr_phi_vec[i])
        GENlep_id = []
        GENlep_Hindex = []
        GENZ_DaughtersId = []
        GENZ_MomId = []
        GENlep_MomId = []
        GENlep_MomMomId = []

        if isMC:
            GENlep_id_vec = self.genworker.GENlep_id
            if len(GENlep_id_vec)>0:
                for i in range(len(GENlep_id_vec)):
                    GENlep_id.append(GENlep_id_vec[i])
            GENlep_Hindex_vec = self.genworker.GENlep_Hindex
            if len(GENlep_Hindex_vec)>0:
                for i in range(len(GENlep_Hindex_vec)):
                    GENlep_Hindex.append(GENlep_Hindex_vec[i])
            GENZ_DaughtersId_vec = self.genworker.GENZ_DaughtersId
            if len(GENZ_DaughtersId_vec)>0:
                for i in range(len(GENZ_DaughtersId_vec)):
                    GENZ_DaughtersId.append(GENZ_DaughtersId_vec[i])
            nVECZ = self.genworker.nVECZ
            GENZ_MomId_vec = self.genworker.GENZ_MomId
            if len(GENZ_MomId_vec)>0:
                for i in range(len(GENZ_MomId_vec)):
                    GENZ_MomId.append(GENZ_MomId_vec[i])
            GENlep_MomId_vec = self.genworker.GENlep_MomId
            if len(GENlep_MomId_vec)>0:
                for i in range(len(GENlep_MomId_vec)):
                    GENlep_MomId.append(GENlep_MomId_vec[i])
            GENlep_MomMomId_vec = self.genworker.GENlep_MomMomId
            if len(GENlep_MomMomId_vec)>0:
                for i in range(len(GENlep_MomMomId_vec)):
                    GENlep_MomMomId.append(GENlep_MomMomId_vec[i])

        if self.worker.GetZ1_2l2qOR2l2nu() and (self.channels == "all"  or self.channels == "2l2v" or self.channels == "2l2q"):  #commented out for now
            if self.channels == "2l2q" or self.channels == "all":
                foundZZCandidate_2l2q = self.worker.ZZSelection_2l2q()
                isBoosted2l2q = self.worker.isBoosted2l2q
                if self.DEBUG: print("isBoosted2l2q: ", isBoosted2l2q)
            if self.channels == "2l2v" or self.channels == "all":
                foundZZCandidate_2l2nu = self.worker.ZZSelection_2l2nu()  #commented out for now
        if (self.channels == "all"  or self.channels == "4l"):
            foundZZCandidate_4l = self.worker.ZZSelection_4l()
            passedFullSelection=foundZZCandidate_4l
            if foundZZCandidate_4l:
                self.passZZEvts += 1
            if (foundZZCandidate_4l |passedFiducialSelection ):
                EvtNum += 1
                keepIt = True
            Lepointer = self.worker.Lepointer
            lep_Hindex = []
            lep_Hindex_vec = self.worker.lep_Hindex
            if len(lep_Hindex_vec)>0:
                for i in range(len(lep_Hindex_vec)):
                    lep_Hindex.append(lep_Hindex_vec[i])
            lep_genindex = []
            if isMC:
                lep_genindex_vec = self.worker.lep_genindex
                if len(lep_genindex_vec)>0:
                    for i in range(len(lep_genindex_vec)):
                        lep_genindex.append(lep_genindex_vec[i])
            if self.worker.RecoFourMuEvent: finalState = 1
            if self.worker.RecoFourEEvent: finalState = 2
            if self.worker.RecoTwoETwoMuEvent: finalState = 3
            if self.worker.RecoTwoMuTwoEEvent: finalState = 4
            if self.worker.flag4e: mass4e = mass4l
            if self.worker.flag2e2mu: mass2e2mu = mass4l
            if self.worker.flag4mu: mass4mu = mass4l

        if self.worker.GetZ1_emuCR() and (self.channels == "all"  or self.channels == "2l2v"):
            foundZZCandidate_2l2nu_emuCR = self.worker.ZZSelection_2l2nu()

        if self.DEBUG:
            print("Found ZZ candidate (4l, 2l2q, 2l2nu): ({}, {}, {})".format(foundZZCandidate_4l, foundZZCandidate_2l2q, foundZZCandidate_2l2nu))

        njets_pt30_eta4p7 = self.worker.njets_pt30_eta4p7
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
        HZZ2l2qNu_cutOppositeChargeFlag = self.worker.HZZ2l2qNu_cutOppositeChargeFlag
        HZZ2l2qNu_nJets = self.worker.HZZ2l2qNu_nJets
        HZZ2l2qNu_nTightBtagJets = self.worker.HZZ2l2qNu_nTightBtagJets
        HZZ2l2qNu_nMediumBtagJets = self.worker.HZZ2l2qNu_nMediumBtagJets
        HZZ2l2qNu_nLooseBtagJets = self.worker.HZZ2l2qNu_nLooseBtagJets

        if (foundZZCandidate_4l or foundZZCandidate_2l2q or foundZZCandidate_2l2nu or foundZZCandidate_2l2nu_emuCR):
            if self.DEBUG: print("Found ZZ candidate (4l, 2l2q, 2l2nu): ({}, {}, {})".format(foundZZCandidate_4l, foundZZCandidate_2l2q, foundZZCandidate_2l2nu))
            pTL1 = self.worker.pTL1
            etaL1 = self.worker.etaL1
            phiL1 = self.worker.phiL1
            massL1 = self.worker.massL1
            pTL2 = self.worker.pTL2
            etaL2 = self.worker.etaL2
            phiL2 = self.worker.phiL2
            massL2 = self.worker.massL2

            if pTL2>pTL1:
                pTL1, pTl2 = pTL2, pTL1
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

        if (foundZZCandidate_2l2q):
            keepIt = True
            foundZZCandidate_2l2q = True
            self.passZZ2l2qEvts += 1
            self.CutFlowTable.Fill(4)

            massZ2_2j = self.worker.Z2_2j.M()
            phiZ2_2j = self.worker.Z2_2j.Phi()
            etaZ2_2j = self.worker.Z2_2j.Eta()
            pTZ2_2j = self.worker.Z2_2j.Pt()
            EneZ2_2j = self.worker.Z2_2j.E()

        if (foundZZCandidate_2l2nu):
            keepIt = True
            foundZZCandidate_2l2nu = True
            self.passZZ2l2nuEvts += 1
            self.CutFlowTable.Fill(5)

            #     FatJet_PNZvsQCD = self.worker.FatJet_PNZvsQCD
            #     self.out.fillBranch("FatJet_PNZvsQCD",FatJet_PNZvsQCD)

        if (foundZZCandidate_2l2nu_emuCR):
            keepIt = True
            passZZ2l2nu_emuCR_Selection = True
            self.passZZ2l2nu_emuCR_Evts += 1
            self.CutFlowTable.Fill(6)

        if (foundZZCandidate_2l2nu or foundZZCandidate_2l2nu_emuCR):
            phiZ2_met = self.worker.Z2_met.Phi()
            pTZ2_met = self.worker.Z2_met.Pt()
            EneZ2_met = self.worker.Z2_met.E()
            MT_2l2nu = self.worker.ZZ_metsystem.Mt()

            HZZ2l2nu_ZZmT = self.worker.ZZ_metsystem.Mt()
            HZZ2l2nu_ZZpT = self.worker.ZZ_metsystem.Pt()

            HZZ2l2nu_ZZmT = self.worker.ZZ_metsystem.Mt()
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
        if (foundZZCandidate_4l):
            keepIt = True
            self.passZZ4lEvts += 1
            self.CutFlowTable.Fill(3)
            foundZZCandidate_4l = True
            print("Inside 4l loop: ",foundZZCandidate_4l)
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
            pTj1 = self.worker.pTj1
            etaj1 = self.worker.etaj1
            phij1 = self.worker.phij1
            mj1 = self.worker.mj1
            pTj2 = self.worker.pTj2
            etaj2 = self.worker.etaj2
            phij2 = self.worker.phij2
            mj2 = self.worker.mj2

            if pTL4>pTL3:
                pTL3, pTL4 = pTL4, pTL3
                etaL3, etaL4 = etaL4, etaL3
                phiL3, phiL4 = phiL4, phiL3
                massL3, massL4 = massL4, massL3

            pT4l = self.worker.ZZsystem.Pt()
            eta4l = self.worker.ZZsystem.Eta()
            phi4l = self.worker.ZZsystem.Phi()
            mass4l = self.worker.ZZsystem.M()
            rapidity4l = self.worker.ZZsystem.Rapidity()
            njets_pt30_eta4p7 = self.worker.njets_pt30_eta4p7
            if self.worker.flag4e:
                mass4e = mass4l
            if self.worker.flag2e2mu:
                mass2e2mu = mass4l
            if self.worker.flag4mu:
                mass4mu = mass4l

            if (self.worker.isFSR==False & passedFullSelection):
                pT4l = self.worker.ZZsystemnofsr.Pt()
                eta4l = self.worker.ZZsystemnofsr.Eta()
                phi4l = self.worker.ZZsystemnofsr.Phi()
                mass4l = self.worker.ZZsystemnofsr.M()
                rapidity4l = self.worker.ZZsystemnofsr.Rapidity()

                Weight = event.genWeight * pileupWeight * dataMCWeight_new * prefiringWeight


        if self.DEBUG:
            print("(found candidates: 2l2q, 2l2nu, 4l): ({:1}, {:1}, {:1}), pTL1: {:7.3f}, pTL2: {:7.3f}, pTZ1: {:7.3f}, pTZ2: {:7.3f}, pTZ2_2j: {:7.3f}, pTZ2_met: {:7.3f}".format(foundZZCandidate_2l2q, foundZZCandidate_2l2nu, foundZZCandidate_4l, pTL1, pTL2, pTZ1, pTZ2, pTZ2_2j, pTZ2_met))
            if (foundZZCandidate_2l2q):
                print("==> pTL1: {}, \t pTL2: {}".format(pTL1, pTL2))

        # Fill the branches with the Trigger information for each channel
        for TriggerChannel in self.cfg['TriggerChannels']:
            self.out.fillBranch(TriggerChannel, TriggerMap[TriggerChannel])
        self.out.fillBranch("mass4l",mass4l)
        self.out.fillBranch("GENmass4l",GENmass4l)
        self.out.fillBranch("mass4e",mass4e)
        self.out.fillBranch("mass2e2mu",mass2e2mu)
        self.out.fillBranch("mass4mu",mass4mu)
        self.out.fillBranch("pT4l",pT4l)
        self.out.fillBranch("GENpT4l",GENpT4l)
        self.out.fillBranch("rapidity4l",rapidity4l)
        self.out.fillBranch("GENrapidity4l",GENrapidity4l)
        self.out.fillBranch("njets_pt30_eta4p7",njets_pt30_eta4p7)
        self.out.fillBranch("finalState",finalState)
        self.out.fillBranch("GENnjets_pt30_eta4p7",GENnjets_pt30_eta4p7)
        self.out.fillBranch("eta4l",eta4l)
        self.out.fillBranch("phi4l",phi4l)
        self.out.fillBranch("massZ1",massZ1)
        self.out.fillBranch("pTZ1",pTZ1)
        self.out.fillBranch("etaZ1",etaZ1)
        self.out.fillBranch("phiZ1",phiZ1)
        self.out.fillBranch("massZ2",massZ2)
        self.out.fillBranch("pTZ2",pTZ2)
        self.out.fillBranch("etaZ2",etaZ2)
        self.out.fillBranch("phiZ2",phiZ2)
        self.out.fillBranch("D_CP",D_CP)
        self.out.fillBranch("D_0m",D_0m)
        self.out.fillBranch("D_0hp",D_0hp)
        self.out.fillBranch("D_int",D_int)
        self.out.fillBranch("D_L1",D_L1)
        self.out.fillBranch("D_L1Zg",D_L1Zg)
        self.out.fillBranch("passedTrig",  passedTrig)
        self.out.fillBranch("passedFullSelection",  passedFullSelection)
        self.out.fillBranch("passedZ4lSelection", passedZ4lSelection)
        self.out.fillBranch("passedZ4lZ1LSelection",  passedZ4lZ1LSelection)
        self.out.fillBranch("passedZ4lZXCRSelection",  passedZ4lZXCRSelection)
        self.out.fillBranch("passedZXCRSelection",  passedZXCRSelection)
        self.out.fillBranch("passedFiducialSelection",  passedFiducialSelection)
        self.out.fillBranch("EvtNum",EvtNum)
        self.out.fillBranch("massL1",massL1)
        self.out.fillBranch("pTL1",pTL1)
        self.out.fillBranch("etaL1",etaL1)
        self.out.fillBranch("phiL1",phiL1)
        self.out.fillBranch("massL2",massL2)
        self.out.fillBranch("pTL2",pTL2)
        self.out.fillBranch("etaL2",etaL2)
        self.out.fillBranch("phiL2",phiL2)
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
        self.out.fillBranch("pileupWeight",pileupWeight)
        self.out.fillBranch("dataMCWeight_new",dataMCWeight_new)
        self.out.fillBranch("prefiringWeight",prefiringWeight)
        self.out.fillBranch("Weight",Weight)
        # self.out.fillBranch("nElectron_Fsr", len(electrons))
        # self.out.fillBranch("nMuon_Fsr", len(muons))

        self.out.fillBranch("GENlep_id",GENlep_id)
        self.out.fillBranch("GENlep_Hindex",GENlep_Hindex)
        self.out.fillBranch("GENZ_DaughtersId",GENZ_DaughtersId)
        self.out.fillBranch("GENZ_MomId",GENZ_MomId)
        self.out.fillBranch("GENlep_MomId",GENlep_MomId)
        self.out.fillBranch("GENlep_MomMomId",GENlep_MomMomId)
        self.out.fillBranch("Electron_Fsr_pt",Electron_Fsr_pt)
        self.out.fillBranch("Electron_Fsr_eta",Electron_Fsr_eta)
        self.out.fillBranch("Electron_Fsr_phi",Electron_Fsr_phi)

        self.out.fillBranch("lep_Hindex",lep_Hindex)
        self.out.fillBranch("lep_genindex",lep_genindex)
        self.out.fillBranch("Muon_Fsr_pt",Muon_Fsr_pt)
        self.out.fillBranch("Muon_Fsr_eta",Muon_Fsr_eta)
        self.out.fillBranch("Muon_Fsr_phi",Muon_Fsr_phi)


        self.out.fillBranch("phiZ2_met",phiZ2_met)
        self.out.fillBranch("pTZ2_met",pTZ2_met)
        self.out.fillBranch("EneZ2_met",EneZ2_met)
        self.out.fillBranch("MT_2l2nu",MT_2l2nu)
        self.out.fillBranch("HZZ2l2nu_ZZmT", HZZ2l2nu_ZZmT)
        self.out.fillBranch("HZZ2l2nu_ZZpT", HZZ2l2nu_ZZpT)
        self.out.fillBranch("HZZ2l2nu_minDPhi_METAK4", HZZ2l2nu_minDPhi_METAK4)
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

        return keepIt


# define modules using the syntax 'name = lambda : constructor' to avoid
# having them loaded when not needed

#H4LCppModule() = lambda: HZZAnalysisCppProducer(year)
