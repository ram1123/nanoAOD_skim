from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import ROOT
import yaml
import os
from modules.Helper import *
from modules.METFilters import passFilters

ROOT.PyConfig.IgnoreCommandLineOptions = True


class HZZAnalysisCppProducer(Module):
    def __init__(self, year, cfgFile, isMC, isFSR, channels, DEBUG=False):
        self.loadLibraries()
        self.year = int(year)
        self.isMC = isMC
        self.channels = channels # choices=["all", "4l", "2l2q", "2l2v"],
        self.DEBUG = DEBUG
        self.cfgFile = cfgFile
        self.cfg = self._load_config(cfgFile)
        self.worker = ROOT.H4LTools(self.year, self.isMC, self.DEBUG)
        self.genworker = ROOT.GenAnalysis(self.DEBUG)
        self._initialize_worker()
        self.worker.isFSR = isFSR
        self._initialize_cutflow_table()
        self._initialize_counters()

        if self.year == 2022:
            self.PUweight_list = self._get_pu_weights()

    def _initialize_cutflow_table(self):
        """Initialize the cut flow table."""
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

        total_cuts = 7 + len(self.dynamicCuts_4l) + len(self.dynamicCuts_2l2q) + len(self.dynamicCuts_2l2nu) + len(self.dynamicCuts_2l2nu_emu_CR)
        print("total_cuts: ", total_cuts)
        self.CutFlowTable =  ROOT.TH1F('cutFlow','cutFlow',total_cuts, 0, total_cuts)
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

    def _get_pu_weights(self):
        """Retrieve pileup weights for the given year."""
        if self.year != 2022:
            return []

        PUinput_file = ROOT.TFile.Open(self.cfg["outputdataNPV"])
        PUinput_hist = PUinput_file.Get(self.cfg["PUweightHistoName"])
        PUweight_list = [PUinput_hist.GetBinContent(i) for i in range(1, PUinput_hist.GetNbinsX() + 1)]
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
            if os.path.exists(fullPath):
                ROOT.gSystem.Load(fullPath)
            else:
                print(f"ERROR: Library {lib} not found at {fullPath}")
                sys.exit(1)

        # Load the yaml-cpp library
        yaml_cpp_lib_path = os.path.join(yaml_cpp_path, "build")
        yaml_cpp_lib = os.path.join(yaml_cpp_lib_path, "libyaml-cpp.so")
        if os.path.exists(yaml_cpp_lib):
            ROOT.gSystem.Load(yaml_cpp_lib)
        else:
            print(f"ERROR: YAML-CPP library not found at {yaml_cpp_lib}")
            sys.exit(1)

        # Load the C++ modules
        h4l_tools_path = os.path.join(base_path, "src/H4LTools.cc")
        gen_analysis_path = os.path.join(base_path, "src/GenAnalysis.cc")

        if "/H4LTools_cc.so" not in ROOT.gSystem.GetLibraries():
            if os.path.exists(h4l_tools_path):
                ROOT.gROOT.ProcessLine(f".L {h4l_tools_path}+O")
            else:
                print(f"Error: H4LTools.cc not found at {h4l_tools_path}")
                sys.exit(1)

        if "/GenAnalysis_cc.so" not in ROOT.gSystem.GetLibraries():
            if os.path.exists(gen_analysis_path):
                ROOT.gROOT.ProcessLine(f".L {gen_analysis_path}+O")
            else:
                print(f"Error: GenAnalysis.cc not found at {gen_analysis_path}")
                sys.exit(1)

    def _load_config(self, cfgFile):
        """Load configuration from YAML file."""
        with open(cfgFile, 'r') as ymlfile:
            return yaml.safe_load(ymlfile)

    def _initialize_worker(self):
        """Initialize the C++ worker with configuration."""
        cfg = self.cfg
        self.worker.InitializeElecut(*self._get_nested_values(cfg['Electron'], [
            'pTcut', 'Etacut', 'Sip3dcut', 'Loosedxycut', 'Loosedzcut', 'Isocut',
            ['BDTWP', 'LowEta', 'LowPT'], ['BDTWP', 'MedEta', 'LowPT'], ['BDTWP', 'HighEta', 'LowPT'],
            ['BDTWP', 'LowEta', 'HighPT'], ['BDTWP', 'MedEta', 'HighPT'], ['BDTWP', 'HighEta', 'HighPT']
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

        self.worker.InitializeEvtCut(*self._get_nested_values(cfg, [
            'MZ1cut', 'MZZcut', ['Higgscut', 'down'], ['Higgscut', 'up'], 'Zmass',
            ['MZcut', 'down'], ['MZcut', 'up'], ['Jet', 'deepJet_btag', 'Loose'],
            ['Jet', 'deepJet_btag', 'Medium'], ['Jet', 'deepJet_btag', 'Tight']
        ]))

        if self.channels in ("all", "2l2q"):
            self.worker.InitializeHZZ2l2qCut(*self._get_nested_values(cfg['HZZ2l2q'], [
                'Leading_Lep_pT', 'SubLeading_Lep_pT', 'Lep_eta', ['MZLepcut', 'down'], ['MZLepcut', 'up']
            ]))

        if self.channels in ("all", "2l2v"):
            self.worker.InitializeHZZ2l2nuCut(*self._get_nested_values(cfg['HZZ2l2nu'], [
                'Leading_Lep_pT', 'SubLeading_Lep_pT', 'Lep_eta', 'Pt_ll', 'M_ll_Window', 'dPhi_jetMET',
                ['MZLepcut', 'down'], ['MZLepcut', 'up']
            ]))


    def _get_nested_values(self, dictionary, keys):
        """Retrieve nested configuration values."""
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
        """Initialize counters for various event types."""
        # self.counters = {
        #     "allEvents": 0,
        #     "triggerPassed": 0,
        #     "metFiltersPassed": 0,
        #     "ZZ4lPassed": 0,
        #     "ZZ2l2qPassed": 0,
        #     "ZZ2l2nuPassed": 0,
        # }
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
        # # print("Processed Events: ", self.counters)
        # print("PassTrig: "+str(self.passtrigEvts)+" Events")
        # print("Pass4eCut: "+str(self.worker.cut4e)+" Events")
        # print("Pass4eGhostRemoval: "+str(self.worker.cutghost4e)+" Events")
        # print("Pass4eLepPtCut: "+str(self.worker.cutLepPt4e)+" Events")
        # print("Pass4eQCDSupress: "+str(self.worker.cutQCD4e)+" Events")
        # print("PassmZ1mZ2Cut_4e: "+str(self.worker.cutZZ4e)+" Events")
        # print("Passm4l_105_160_Cut_4e: "+str(self.worker.cutm4l4e)+" Events")
        # print("Pass4muCut: "+str(self.worker.cut4mu)+" Events")
        # print("Pass4muGhostRemoval: "+str(self.worker.cutghost4mu)+" Events")
        # print("Pass4muLepPtCut: "+str(self.worker.cutLepPt4mu)+" Events")
        # print("Pass4muQCDSupress: "+str(self.worker.cutQCD4mu)+" Events")
        # print("PassmZ1mZ2Cut_4mu: "+str(self.worker.cutZZ4mu)+" Events")
        # print("Passm4l_105_160_Cut_4mu: "+str(self.worker.cutm4l4mu)+" Events")
        # print("Pass2e2muCut: "+str(self.worker.cut2e2mu)+" Events")
        # print("Pass2e2muGhostRemoval: "+str(self.worker.cutghost2e2mu)+" Events")
        # print("Pass2e2muLepPtCut: "+str(self.worker.cutLepPt2e2mu)+" Events")
        # print("Pass2e2muQCDSupress: "+str(self.worker.cutQCD2e2mu)+" Events")
        # print("PassmZ1mZ2Cut_2e2mu: "+str(self.worker.cutZZ2e2mu)+" Events")
        # print("Passm4l_105_160_Cut_2e2mu: "+str(self.worker.cutm4l2e2mu)+" Events")
        # # print("PassZZSelection: "+str(self.passZZEvts)+" Events")
        # if self.isMC and self.year == 2022:
        #     print("PassGEN4eCut: "+str(self.genworker.nGEN4e)+" Events")
        #     print("PassGEN4eZ1Cut: "+str(self.genworker.nGEN4epassZ1)+" Events")
        #     print("PassGEN4efidCut: "+str(self.genworker.nGEN4epassFid)+" Events")
        #     print("PassGEN2e2muCut: "+str(self.genworker.nGEN2e2mu)+" Events")
        #     print("PassGEN2e2muZ1Cut: "+str(self.genworker.nGEN2e2mupassZ1)+" Events")
        #     print("PassGEN2e2mufidCut: "+str(self.genworker.nGEN2e2mupassFid)+" Events")
        #     print("PassGEN4muCut: "+str(self.genworker.nGEN4mu)+" Events")
        #     print("PassGEN4muZ1Cut: "+str(self.genworker.nGEN4mupassZ1)+" Events")
        #     print("PassGEN4mufidCut: "+str(self.genworker.nGEN4mupassFid)+" Events")
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree)  # initReaders must be called in beginFile
        self.out = wrappedOutputTree

        # Boolean branches for Trigger channels
        for TriggerChannel in self.cfg['TriggerChannels']:
            self.out.branch(TriggerChannel, "O", title=f"Trigger channel: {TriggerChannel}")

        # Boolean branches for channel candidates
        self.out.branch("foundZZCandidate_4l", "O", title="ZZ Candidate found in 4l channel")
        self.out.branch("foundZZCandidate_2l2q", "O", title="ZZ Candidate found in 2l2q channel")
        self.out.branch("foundZZCandidate_2l2nu", "O", title="ZZ Candidate found in 2l2nu channel")
        self.out.branch("passZZ2l2nu_emuCR_Selection", "O", title="Pass 2l2nu emu control region selection")
        self.out.branch("isBoosted2l2q", "O", title="Boosted topology in 2l2q channel")
        self.out.branch("HZZ2l2nu_ifVBF", "O", title="VBF topology in 2l2nu channel")
        self.out.branch("HZZ2l2qNu_isELE", "O", title="Electron topology in 2l2q channel")
        self.out.branch("HZZ2l2qNu_cutOppositeChargeFlag", "O", title="Opposite charge cut passed for 2l2qNu")
        self.out.branch("HZZ2l2nu_isEMuCR", "O", title="EMU control region in 2l2nu channel")

        # Branches for lepton-related variables
        for i in range(1, 5):
            self.out.branch(f"massL{i}", "F", title=f"Mass of lepton {i}")
            self.out.branch(f"pTL{i}", "F", title=f"Transverse momentum (pT) of lepton {i}")
            self.out.branch(f"etaL{i}", "F", title=f"Pseudo-rapidity (eta) of lepton {i}")
            self.out.branch(f"phiL{i}", "F", title=f"Azimuthal angle (phi) of lepton {i}")

        # Branches for 4l channel: ZZ kinematics
        self.out.branch("mass4l", "F", title="Mass of ZZ system in 4l channel")
        self.out.branch("GENmass4l", "F", title="Generated mass of ZZ system in 4l channel")
        self.out.branch("mass4e", "F", title="Mass of 4 electrons")
        self.out.branch("mass4mu", "F", title="Mass of 4 muons")
        self.out.branch("mass2e2mu", "F", title="Mass of 2 electrons and 2 muons")
        self.out.branch("pT4l", "F", title="Transverse momentum (pT) of ZZ system")
        self.out.branch("GENpT4l", "F", title="Generated transverse momentum (pT) of ZZ system")
        self.out.branch("rapidity4l", "F", title="Rapidity of ZZ system")
        self.out.branch("njets_pt30_eta4p7", "I", title="Number of jets with pT > 30 GeV and |eta| < 4.7")
        self.out.branch("finalState", "I", title="Final state identifier for 4l channel")
        self.out.branch("GENnjets_pt30_eta4p7", "I", title="Generated number of jets with pT > 30 GeV and |eta| < 4.7")
        self.out.branch("GENrapidity4l", "F", title="Generated rapidity of ZZ system")
        self.out.branch("eta4l", "F", title="Pseudo-rapidity (eta) of ZZ system")
        self.out.branch("phi4l", "F", title="Azimuthal angle (phi) of ZZ system")

        # Branches for 4l channel: MELA discriminants
        mela_discriminants = ["D_CP", "D_0m", "D_0hp", "D_int", "D_L1", "D_L1Zg"]
        for discriminant in mela_discriminants:
            self.out.branch(discriminant, "F", title=f"MELA discriminant: {discriminant}")

        # Branches for jets in the 4l channel
        for i in range(1, 3):
            self.out.branch(f"mj{i}", "F", title=f"Mass of jet {i}")
            self.out.branch(f"pTj{i}", "F", title=f"Transverse momentum (pT) of jet {i}")
            self.out.branch(f"etaj{i}", "F", title=f"Pseudo-rapidity (eta) of jet {i}")
            self.out.branch(f"phij{i}", "F", title=f"Azimuthal angle (phi) of jet {i}")

        # Branches common to 4l, 2l2q, and 2l2nu channels
        for z in range(1, 3):
            self.out.branch(f"massZ{z}", "F", title=f"Mass of Z boson {z}")
            self.out.branch(f"pTZ{z}", "F", title=f"Transverse momentum (pT) of Z boson {z}")
            self.out.branch(f"etaZ{z}", "F", title=f"Pseudo-rapidity (eta) of Z boson {z}")
            self.out.branch(f"phiZ{z}", "F", title=f"Azimuthal angle (phi) of Z boson {z}")

        # Additional general branches
        self.out.branch("EvtNum", "I", title="Event number")
        self.out.branch("Weight", "F", title="Event weight")
        self.out.branch("pileupWeight", "F", title="Pileup reweighting factor")
        self.out.branch("dataMCWeight_new", "F", title="Data/MC weight")
        self.out.branch("prefiringWeight", "F", title="Prefiring weight correction")
        self.out.branch("passedTrig", "O", title="Event passed trigger")
        self.out.branch("passedFullSelection", "O", title="Event passed full selection criteria")
        self.out.branch("passedZ4lSelection", "O", title="Event passed 4l ZZ selection")
        self.out.branch("passedZ4lZ1LSelection", "O", title="Event passed Z1+lepton selection in 4l channel")
        self.out.branch("passedZ4lZXCRSelection", "O", title="Event passed Z+X control region selection")
        self.out.branch("passedZXCRSelection", "O", title="Event passed ZX control region selection")
        self.out.branch("passedFiducialSelection", "O", title="Event passed fiducial selection")

        # GEN-level branches
        GENHlepNum = 4
        GENZNum = 2
        self.out.branch("GENlep_MomId", "I", lenVar="nGenPart", title="GEN lepton mother ID")
        self.out.branch("GENlep_MomMomId", "I", lenVar="nGenPart", title="GEN lepton grandmother ID")
        self.out.branch("GENZ_MomId", "I", lenVar="nVECZ", title="GEN Z boson mother ID")
        self.out.branch("GENZ_DaughtersId", "I", lenVar="GENZNum", title="GEN Z boson daughter IDs")
        self.out.branch("GENlep_Hindex", "I", lenVar="GENHlepNum", title="GEN lepton Higgs indices")
        self.out.branch("lep_Hindex", "I", lenVar="GENHlepNum", title="Reco lepton Higgs indices")
        self.out.branch("GENlep_id", "I", lenVar="nGENLeptons", title="GEN lepton IDs")
        self.out.branch("lep_genindex", "I", lenVar="Lepointer", title="Reco lepton GEN indices")

        # FSR branches for leptons
        for particle in ["Electron", "Muon"]:
            self.out.branch(f"{particle}_Fsr_pt", "F", lenVar=f"n{particle}_Fsr", title=f"{particle} FSR photon transverse momentum")
            self.out.branch(f"{particle}_Fsr_eta", "F", lenVar=f"n{particle}_Fsr", title=f"{particle} FSR photon pseudo-rapidity")
            self.out.branch(f"{particle}_Fsr_phi", "F", lenVar=f"n{particle}_Fsr", title=f"{particle} FSR photon azimuthal angle")

        # Branches for 2l2q channel
        for var in ["massZ2_2j", "phiZ2_2j", "etaZ2_2j", "pTZ2_2j", "EneZ2_2j", "phiZ2_met", "pTZ2_met", "EneZ2_met", "MT_2l2nu"]:
            self.out.branch(var, "F", title=f"Variable {var} for 2l2q/2l2nu channel")


        # VBF jets and dijet kinematics for 2l2nu channel
        for var in [ "HZZ2l2qNu_nJets", "HZZ2l2qNu_nTightBtagJets", "HZZ2l2qNu_nMediumBtagJets", "HZZ2l2qNu_nLooseBtagJets"]:
            self.out.branch(var, "I", title=f"Variable {var} for VBF jets")

        for jet in ["VBFjet1", "VBFjet2", "VBFdijet"]:
            for attr in ["pT", "eta", "phi", "mass", "E"]:
                self.out.branch(f"HZZ2l2nu_{jet}_{attr}", "F", title=f"{jet} {attr} for 2l2nu channel")

        for var in ["HZZ2l2nu_VBFIndexJet1", "HZZ2l2nu_VBFIndexJet2"]:
            self.out.branch(var, "I", title=f"Variable {var} for VBF jets")

        for var in ["HZZ2l2nu_ZZmT", "HZZ2l2nu_ZZpT", "HZZ2l2nu_minDPhi_METAK4", "HZZ2l2nu_VBFdEta_jj", "HZZ2l2nu_VBFdPhi_jj", "HZZ2l2nu_VBFdR_jj"]:
            self.out.branch(var, "F", title=f"Variable {var} for VBF jets")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
         print("\n========== Print Cut flow table  ====================\n")
         self.cutFlowCounts = {
             "Total": self.passAllEvts,
             "PassTrig": self.passtrigEvts,
             "PassMETFilters": self.passMETFilters,
             "PassZZSelection": self.passZZ4lEvts,


             "PassZZ2l2qSelection": self.passZZ2l2qEvts,
             "PassZZ2l2nuSelection": self.passZZ2l2nuEvts,
             "PassZZ2l2nu_emuCR_Selection": self.passZZ2l2nu_emuCR_Evts
         }

         # Cut flow data for 4l, 2l2q, 2l2nu channels for json output
         cutFlowData = {
             "General": self.cutFlowCounts,
             "4l_Channel": {},
             "2l2q_Channel": {},
             "2l2nu_Channel": {},
             "2l2nu_Channel_emu_CR": {}
         }

         for key, value in self.cutFlowCounts.items():
             print("{:27}:{:7} {}".format(key, str(value), " Events"))

         print("\n4l channel cut flow table:")
         for idx, cut in enumerate(self.dynamicCuts_4l):
             print("{:2} {:27}:{:7} {}".format(8 + idx, cut, str(getattr(self.worker, cut)), " Events"))
             cutFlowData["4l_Channel"][cut] = getattr(self.worker, cut, 'N/A')
             self.CutFlowTable.SetBinContent(8 + idx, getattr(self.worker, cut, 'N/A'))

         print("\n2l2q channel cut flow table:")
         for idx, cut in enumerate(self.dynamicCuts_2l2q):
             print("{:2} {:27}:{:7} {}".format(8 + len(self.dynamicCuts_4l) + idx, cut, str(getattr(self.worker, cut)), " Events"))
             cutFlowData["2l2q_Channel"][cut] = getattr(self.worker, cut, 'N/A')
             self.CutFlowTable.SetBinContent(8 + len(self.dynamicCuts_4l) + idx, getattr(self.worker, cut, 'N/A'))

         print("\n2l2nu channel cut flow table:")
         for idx, cut in enumerate(self.dynamicCuts_2l2nu):
             print("{:2} {:27}:{:7} {}".format(8 + len(self.dynamicCuts_4l) + len(self.dynamicCuts_2l2q) + idx, cut, str(getattr(self.worker, cut)), " Events"))
             cutFlowData["2l2nu_Channel"][cut] = getattr(self.worker, cut, 'N/A')
             self.CutFlowTable.SetBinContent(8 + len(self.dynamicCuts_4l) + len(self.dynamicCuts_2l2q) + idx, getattr(self.worker, cut, 'N/A'))

         print("\n2l2nu channel emu control region cut flow table:")
         for idx, cut in enumerate(self.dynamicCuts_2l2nu_emu_CR):
             print("{:2} {:27}:{:7} {}".format(8 + len(self.dynamicCuts_4l) + len(self.dynamicCuts_2l2q) + len(self.dynamicCuts_2l2nu) + idx, cut, str(getattr(self.worker, cut)), " Events"))
             cutFlowData["2l2nu_Channel_emu_CR"][cut] = getattr(self.worker, cut, 'N/A')
             self.CutFlowTable.SetBinContent(8 + len(self.dynamicCuts_4l) + len(self.dynamicCuts_2l2q) + len(self.dynamicCuts_2l2nu) + idx, getattr(self.worker, cut, 'N/A'))

         print("\n========== END: Print Cut flow table  ====================\n")

         outputFile.cd()
         self.CutFlowTable.Write()

    # this function gets the pointers to Value and ArrayReaders and sets
    # them in the C++ worker class
    def initReaders(self, tree):
        # self._ttreereaderversion must be set AFTER all calls to
        # tree.valueReader or tree.arrayReader
        self._ttreereaderversion = tree._ttreereaderversion


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail,
        go to next event)"""
        keepIt = False
        self.worker.Initialize()
        self.passAllEvts += 1
        self.CutFlowTable.Fill(0)

        # Initialize variables for ZZ candidates
        pTZ1 = -99
        etaZ1 = -99
        phiZ1 = -99
        massZ1 = 0
        pTZ2 = -99
        etaZ2 = -99
        phiZ2 = -99
        massZ2 = -99

        # Initialize variables for lepton properties
        pTL1 = -99
        etaL1 = -99
        phiL1 = -99
        massL1 = -99
        pTL2 = -99
        etaL2 = -99
        phiL2 = -99
        massL2 = -99
        pTL3 = -99
        etaL3 = -99
        phiL3 = -99
        massL3 = -99
        pTL4 = -99
        etaL4 = -99
        phiL4 = -99
        massL4 = -99

        # Initialize variables for ZZ system properties
        mass4l = -99
        GENmass4l = -99
        mass4e = -99
        mass4mu = -99
        mass2e2mu = -99
        pT4l = -99
        GENpT4l = -99
        rapidity4l = -99
        GENrapidity4l = -99
        eta4l = -99
        phi4l = -99
        njets_pt30_eta4p7 = -1
        GENnjets_pt30_eta4p7 = -1

        # Initialize MELA discriminants
        D_CP = -99
        D_0m = -99
        D_0hp = -99
        D_int = -99
        D_L1 = -99
        D_L1Zg = -99

        # Initialize jet properties for 4l channel
        pTj1 = -99
        etaj1 = -99
        phij1 = -99
        mj1 = -99
        pTj2 = -99
        etaj2 = -99
        phij2 = -99
        mj2 = -99

        # Initialize 2l2q-specific variables
        massZ2_2j = -99
        phiZ2_2j = -99
        etaZ2_2j = -99
        pTZ2_2j = -99
        EneZ2_2j = -99

        # Initialize 2l2nu-specific variables
        phiZ2_met = -99
        pTZ2_met = -99
        EneZ2_met = -99
        MT_2l2nu = -99
        HZZ2l2nu_ZZmT = -99
        HZZ2l2nu_ZZpT = -99
        HZZ2l2nu_minDPhi_METAK4 = -99

        # Initialize VBF jet and dijet properties
        HZZ2l2nu_VBFIndexJet1 = -1
        HZZ2l2nu_VBFIndexJet2 = -1
        HZZ2l2nu_VBFjet1_pT = -99
        HZZ2l2nu_VBFjet1_eta = -99
        HZZ2l2nu_VBFjet1_phi = -99
        HZZ2l2nu_VBFjet1_mass = -99
        HZZ2l2nu_VBFjet2_pT = -99
        HZZ2l2nu_VBFjet2_eta = -99
        HZZ2l2nu_VBFjet2_phi = -99
        HZZ2l2nu_VBFjet2_mass = -99
        HZZ2l2nu_VBFdijet_mass = -99
        HZZ2l2nu_VBFdijet_pT = -99
        HZZ2l2nu_VBFdijet_E = -99
        HZZ2l2nu_VBFdEta_jj = -99
        HZZ2l2nu_VBFdPhi_jj = -99
        HZZ2l2nu_VBFdR_jj = -99

        # Initialize GEN-level variables
        GENlep_MomId = []
        GENlep_MomMomId = []
        GENlep_id = []
        GENlep_Hindex = []
        lep_Hindex = []
        lep_genindex = []
        GENZ_DaughtersId = []
        GENZ_MomId = []
        Electron_Fsr_pt = []
        Electron_Fsr_eta = []
        Electron_Fsr_phi = []
        Muon_Fsr_pt = []
        Muon_Fsr_eta = []
        Muon_Fsr_phi = []

        # Initialize event-level variables
        EvtNum = 0
        Weight = 1
        pileupWeight = 1
        dataMCWeight_new = 1
        prefiringWeight = 1
        finalState = -1

        # Initialize flags
        foundZZCandidate_4l = False
        foundZZCandidate_2l2q = False
        foundZZCandidate_2l2nu = False
        foundZZCandidate_2l2nu_emuCR = False   # for 2l2nu emu control region
        passZZ2l2nu_emuCR_Selection = False
        isBoosted2l2q = False
        HZZ2l2nu_ifVBF = False
        HZZ2l2qNu_isELE = False
        HZZ2l2nu_isEMuCR = False
        HZZ2l2qNu_cutOppositeChargeFlag = False
        passedTrig = False
        passedFullSelection = False
        passedZ4lSelection = False
        passedZ4lZ1LSelection = False
        passedZ4lZXCRSelection = False
        passedZXCRSelection = False
        passedFiducialSelection = False

        # Other variables
        GENHlepNum = 4
        GENZNum = 2

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

        # Pass MET filters
        if not passFilters(event, self.year):
            return keepIt

        self.passMETFilters += 1
        self.CutFlowTable.Fill(2)

        # Pileup weight
        pileupWeight = self.PUweight_list[event.Pileup_nPU] if self.year == 2022 and self.isMC else 1.0

        # Process objects
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        # Photons = Collection(event, "Photon")
        jets = Collection(event, "Jet")
        FatJets = Collection(event, "FatJet")
        met = Object(event, "MET", None)

        if self.isMC and self.year == 2022:
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

            lep_genindex = []
            lep_genindex_vec = self.worker.lep_genindex
            if len(lep_genindex_vec)>0:
                for i in range(len(lep_genindex_vec)):
                    lep_genindex.append(lep_genindex_vec[i])

            GENlep_id = []
            GENlep_Hindex = []
            GENZ_DaughtersId = []
            GENZ_MomId = []
            GENlep_MomId = []
            GENlep_MomMomId = []
            self.genworker.SetGenVariables()
            GENmass4l = self.genworker.GENmass4l
            GENpT4l = self.genworker.GENpT4l
            GENrapidity4l = self.genworker.GENrapidity4l
            GENnjets_pt30_eta4p7 = self.genworker.GENnjets_pt30_eta4p7
            nGENLeptons = self.genworker.nGENLeptons
            passedFiducialSelection = self.genworker.passedFiducialSelection

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

                genWeight = 1 if event.genWeight > 0 else -1
                Weight = genWeight * pileupWeight * dataMCWeight_new * prefiringWeight


        if self.DEBUG:
            print("(found candidates: 2l2q, 2l2nu, 4l): ({:1}, {:1}, {:1}), pTL1: {:7.3f}, pTL2: {:7.3f}, pTZ1: {:7.3f}, pTZ2: {:7.3f}, pTZ2_2j: {:7.3f}, pTZ2_met: {:7.3f}".format(foundZZCandidate_2l2q, foundZZCandidate_2l2nu, foundZZCandidate_4l, pTL1, pTL2, pTZ1, pTZ2, pTZ2_2j, pTZ2_met))
            if (foundZZCandidate_2l2q):
                print("==> pTL1: {}, \t pTL2: {}".format(pTL1, pTL2))

        # Filling the branches with values
        # Fill the branches with the Trigger information for each channel
        for TriggerChannel in self.cfg['TriggerChannels']:
            self.out.fillBranch(TriggerChannel, TriggerMap[TriggerChannel])

        self.out.fillBranch("mass4l", mass4l)
        self.out.fillBranch("GENmass4l", GENmass4l)
        self.out.fillBranch("mass4e", mass4e)
        self.out.fillBranch("mass2e2mu", mass2e2mu)
        self.out.fillBranch("mass4mu", mass4mu)
        self.out.fillBranch("pT4l", pT4l)
        self.out.fillBranch("GENpT4l", GENpT4l)
        self.out.fillBranch("rapidity4l", rapidity4l)
        self.out.fillBranch("GENrapidity4l", GENrapidity4l)
        self.out.fillBranch("njets_pt30_eta4p7", njets_pt30_eta4p7)
        self.out.fillBranch("finalState", finalState)
        self.out.fillBranch("GENnjets_pt30_eta4p7", GENnjets_pt30_eta4p7)
        self.out.fillBranch("eta4l", eta4l)
        self.out.fillBranch("phi4l", phi4l)
        self.out.fillBranch("massZ1", massZ1)
        self.out.fillBranch("pTZ1", pTZ1)
        self.out.fillBranch("etaZ1", etaZ1)
        self.out.fillBranch("phiZ1", phiZ1)
        self.out.fillBranch("massZ2", massZ2)
        self.out.fillBranch("pTZ2", pTZ2)
        self.out.fillBranch("etaZ2", etaZ2)
        self.out.fillBranch("phiZ2", phiZ2)
        self.out.fillBranch("D_CP", D_CP)
        self.out.fillBranch("D_0m", D_0m)
        self.out.fillBranch("D_0hp", D_0hp)
        self.out.fillBranch("D_int", D_int)
        self.out.fillBranch("D_L1", D_L1)
        self.out.fillBranch("D_L1Zg", D_L1Zg)
        self.out.fillBranch("passedTrig", passedTrig)
        self.out.fillBranch("passedFullSelection", passedFullSelection)
        self.out.fillBranch("passedZ4lSelection", passedZ4lSelection)
        self.out.fillBranch("passedZ4lZ1LSelection", passedZ4lZ1LSelection)
        self.out.fillBranch("passedZ4lZXCRSelection", passedZ4lZXCRSelection)
        self.out.fillBranch("passedZXCRSelection", passedZXCRSelection)
        self.out.fillBranch("passedFiducialSelection", passedFiducialSelection)
        self.out.fillBranch("EvtNum", EvtNum)
        self.out.fillBranch("massL1", massL1)
        self.out.fillBranch("pTL1", pTL1)
        self.out.fillBranch("etaL1", etaL1)
        self.out.fillBranch("phiL1", phiL1)
        self.out.fillBranch("massL2", massL2)
        self.out.fillBranch("pTL2", pTL2)
        self.out.fillBranch("etaL2", etaL2)
        self.out.fillBranch("phiL2", phiL2)
        self.out.fillBranch("massL3", massL3)
        self.out.fillBranch("pTL3", pTL3)
        self.out.fillBranch("etaL3", etaL3)
        self.out.fillBranch("phiL3", phiL3)
        self.out.fillBranch("massL4", massL4)
        self.out.fillBranch("pTL4", pTL4)
        self.out.fillBranch("etaL4", etaL4)
        self.out.fillBranch("phiL4", phiL4)
        self.out.fillBranch("mj1", mj1)
        self.out.fillBranch("pTj1", pTj1)
        self.out.fillBranch("etaj1", etaj1)
        self.out.fillBranch("phij1", phij1)
        self.out.fillBranch("mj2", mj2)
        self.out.fillBranch("pTj2", pTj2)
        self.out.fillBranch("etaj2", etaj2)
        self.out.fillBranch("phij2", phij2)
        self.out.fillBranch("pileupWeight", pileupWeight)
        self.out.fillBranch("dataMCWeight_new", dataMCWeight_new)
        self.out.fillBranch("prefiringWeight", prefiringWeight)
        self.out.fillBranch("Weight", Weight)
        self.out.fillBranch("GENlep_id", GENlep_id)
        self.out.fillBranch("GENlep_Hindex", GENlep_Hindex)
        self.out.fillBranch("GENZ_DaughtersId", GENZ_DaughtersId)
        self.out.fillBranch("GENZ_MomId", GENZ_MomId)
        self.out.fillBranch("GENlep_MomId", GENlep_MomId)
        self.out.fillBranch("GENlep_MomMomId", GENlep_MomMomId)
        self.out.fillBranch("Electron_Fsr_pt", Electron_Fsr_pt)
        self.out.fillBranch("Electron_Fsr_eta", Electron_Fsr_eta)
        self.out.fillBranch("Electron_Fsr_phi", Electron_Fsr_phi)
        self.out.fillBranch("lep_Hindex", lep_Hindex)
        self.out.fillBranch("lep_genindex", lep_genindex)
        self.out.fillBranch("Muon_Fsr_pt", Muon_Fsr_pt)
        self.out.fillBranch("Muon_Fsr_eta", Muon_Fsr_eta)
        self.out.fillBranch("Muon_Fsr_phi", Muon_Fsr_phi)
        self.out.fillBranch("phiZ2_met", phiZ2_met)
        self.out.fillBranch("pTZ2_met", pTZ2_met)
        self.out.fillBranch("EneZ2_met", EneZ2_met)
        self.out.fillBranch("MT_2l2nu", MT_2l2nu)
        self.out.fillBranch("HZZ2l2nu_ZZmT", HZZ2l2nu_ZZmT)
        self.out.fillBranch("HZZ2l2nu_ZZpT", HZZ2l2nu_ZZpT)
        self.out.fillBranch("HZZ2l2nu_minDPhi_METAK4", HZZ2l2nu_minDPhi_METAK4)
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
