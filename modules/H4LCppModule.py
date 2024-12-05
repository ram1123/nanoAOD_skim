from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import ROOT
import yaml
import os
from modules.Helper import *
from modules.METFilters import passFilters
from config.branch_config import branch_definitions

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
        self.branch_definitions = branch_definitions
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
        self.passZZEvts = 0
        self.EvtNum = 0
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

        for branch_name, details in self.branch_definitions.items():
            if ("GEN" in branch_name or "lep_genindex" in branch_name) and self.year != 2022:
                # FIXME: need to fix this for other years
                continue
            # if "GEN" in branch_name and self.year == 2022: continue
            if "lenVar" in details:
                if self.DEBUG: print(f"Defining the branch: {branch_name} {details['type']} {details['lenVar']}")
                self.out.branch(branch_name, details["type"], lenVar=f'{details["lenVar"]}', title=details.get("title", ""))
            else:
                self.out.branch(branch_name, details["type"], title=details.get("title", ""))

        if self.DEBUG: print("Branches are defined")

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
        if self.DEBUG: print("\n\n=====================================")
        if self.DEBUG: print("Processing event: ", event.event)
        self.passAllEvts += 1
        if self.DEBUG: print("Event number: ", self.passAllEvts)
        keepIt = False
        self.worker.Initialize()
        if self.isMC: self.worker.SetObjectNumGen(event.nGenPart)
        self.worker.SetObjectNum(event.nElectron,event.nMuon,event.nJet,event.nFsrPhoton)


        self.CutFlowTable.Fill(0)

        self.branch_values = {branch: details["default"] for branch, details in self.branch_definitions.items()}
        if self.DEBUG: print("Branches are initialized")
        if self.DEBUG: print(f"branch values: {self.branch_values}")
        if self.DEBUG: print(f'lep_genindex: {self.branch_values["lep_genindex"]}')
        if self.DEBUG: print(f'GENlep_MomMomId: {self.branch_values["GENlep_MomMomId"]}')

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
        if (self.DEBUG): print("MET filters passed")

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
        if self.DEBUG: print("Objects are collected")

        if self.isMC and self.year == 2022:
            self.genworker.Initialize();
            self.genworker.SetObjectNumGen(event.nGenPart, event.nGenJet)

            if self.DEBUG: print("Processing GEN particles")
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

            self.branch_values["lep_genindex"] = []
            lep_genindex_vec = self.worker.lep_genindex
            if len(lep_genindex_vec)>0:
                print("RAM: len(lep_genindex_vec): ", len(lep_genindex_vec))
                for i in range(len(lep_genindex_vec)):
                    self.branch_values["lep_genindex"].append(lep_genindex_vec[i])

            if self.DEBUG: print("Setting GEN variables")
            self.genworker.SetGenVariables()
            self.branch_values["GENlep_id"] = []
            self.branch_values["GENlep_Hindex"] = []
            self.branch_values["GENZ_DaughtersId"] = []
            self.branch_values["GENZ_MomId"] = []
            self.branch_values["GENlep_MomId"] = []
            self.branch_values["GENlep_MomMomId"] = []
            self.branch_values["GENmass4l"] = self.genworker.GENmass4l
            self.branch_values["GENpT4l"] = self.genworker.GENpT4l
            self.branch_values["GENrapidity4l"] = self.genworker.GENrapidity4l
            self.branch_values["GENnjets_pt30_eta4p7"] = self.genworker.GENnjets_pt30_eta4p7
            self.branch_values["nGENLeptons"] = self.genworker.nGENLeptons
            self.branch_values["passedFiducialSelection"] = self.genworker.passedFiducialSelection
            if self.DEBUG:print("len(self.branch_values[GENlep_MomMomId]): ", len(self.branch_values["GENlep_MomMomId"]))

            if self.DEBUG: print("len(GENlep_id_vec): ", len(self.genworker.GENlep_id))
            GENlep_id_vec = self.genworker.GENlep_id
            if len(GENlep_id_vec)>0:
                for i in range(len(GENlep_id_vec)):
                    self.branch_values["GENlep_id"].append(GENlep_id_vec[i])
            GENlep_Hindex_vec = self.genworker.GENlep_Hindex
            if (self.DEBUG): print("len(GENlep_Hindex_vec): ", len(GENlep_Hindex_vec))
            if len(GENlep_Hindex_vec)>0:
                for i in range(len(GENlep_Hindex_vec)):
                    self.branch_values["GENlep_Hindex"].append(GENlep_Hindex_vec[i])
            GENZ_DaughtersId_vec = self.genworker.GENZ_DaughtersId
            if len(GENZ_DaughtersId_vec)>0:
                for i in range(len(GENZ_DaughtersId_vec)):
                    self.branch_values["GENZ_DaughtersId"].append(GENZ_DaughtersId_vec[i])
            nVECZ = self.genworker.nVECZ
            GENZ_MomId_vec = self.genworker.GENZ_MomId
            if len(GENZ_MomId_vec)>0:
                for i in range(len(GENZ_MomId_vec)):
                    self.branch_values["GENZ_MomId"].append(GENZ_MomId_vec[i])
            GENlep_MomId_vec = self.genworker.GENlep_MomId
            if len(GENlep_MomId_vec)>0:
                for i in range(len(GENlep_MomId_vec)):
                    self.branch_values["GENlep_MomId"].append(GENlep_MomId_vec[i])
            GENlep_MomMomId_vec = self.genworker.GENlep_MomMomId
            if self.DEBUG: print("CHECK SIZE: len(GENlep_MomMomId_vec): ", len(GENlep_MomMomId_vec))
            if len(GENlep_MomMomId_vec)>0:
                for i in range(len(GENlep_MomMomId_vec)):
                    self.branch_values["GENlep_MomMomId"].append(GENlep_MomMomId_vec[i])
            if self.DEBUG: print("GEN variables are set")

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
        # if self.year == 2022:
        self.worker.BatchFsrRecovery_Run3() # It should run for all years, as some of variables are used in the worker. So, if we comment it for other years it will give error

        self.worker.LeptonSelection()
        if self.DEBUG:
            print("Lepton selection done")

        Electron_Fsr_pt_vec = self.worker.ElectronFsrPt()
        Electron_Fsr_eta_vec = self.worker.ElectronFsrEta()
        Electron_Fsr_phi_vec = self.worker.ElectronFsrPhi()
        Muon_Fsr_pt_vec = self.worker.MuonFsrPt()
        Muon_Fsr_eta_vec = self.worker.MuonFsrEta()
        Muon_Fsr_phi_vec = self.worker.MuonFsrPhi()

        self.branch_values["Electron_Fsr_pt"] = []
        self.branch_values["Electron_Fsr_eta"] = []
        self.branch_values["Electron_Fsr_phi"] = []
        self.branch_values["Muon_Fsr_pt"] = []
        self.branch_values["Muon_Fsr_eta"] = []
        self.branch_values["Muon_Fsr_phi"] = []
        self.branch_values["lep_Hindex"] = []

        if len(Electron_Fsr_pt_vec)>0:
            for i in range(len(Electron_Fsr_pt_vec)):
                self.branch_values["Electron_Fsr_pt"].append(Electron_Fsr_pt_vec[i])
                self.branch_values["Electron_Fsr_eta"].append(Electron_Fsr_eta_vec[i])
                self.branch_values["Electron_Fsr_phi"].append(Electron_Fsr_phi_vec[i])
        if len(Muon_Fsr_pt_vec)>0:
            for i in range(len(Muon_Fsr_pt_vec)):
                self.branch_values["Muon_Fsr_pt"].append(Muon_Fsr_pt_vec[i])
                self.branch_values["Muon_Fsr_eta"].append(Muon_Fsr_eta_vec[i])
                self.branch_values["Muon_Fsr_phi"].append(Muon_Fsr_phi_vec[i])

        if self.worker.GetZ1_2l2qOR2l2nu() and (self.channels == "all"  or self.channels == "2l2v" or self.channels == "2l2q"):  #commented out for now
            if self.channels == "2l2q" or self.channels == "all":
                self.branch_values["foundZZCandidate_2l2q"] = self.worker.ZZSelection_2l2q()
                self.branch_values["isBoosted2l2q"] = self.worker.isBoosted2l2q
                if self.DEBUG: print("isBoosted2l2q: ", self.branch_values["isBoosted2l2q"])
            if self.channels == "2l2v" or self.channels == "all":
                self.branch_values["foundZZCandidate_2l2nu"] = self.worker.ZZSelection_2l2nu()  #commented out for now
        if (self.channels == "all"  or self.channels == "4l"):
            self.branch_values["foundZZCandidate_4l"] = self.worker.ZZSelection_4l()
            self.branch_values["passedFullSelection"]= self.branch_values["foundZZCandidate_4l"]
            if  self.branch_values["foundZZCandidate_4l"]:
                self.passZZEvts += 1
            if (self.branch_values["foundZZCandidate_4l"] |self.branch_values["passedFiducialSelection"] ):
                self.EvtNum += 1
                keepIt = True
            Lepointer = self.worker.Lepointer
            lep_Hindex_vec = self.worker.lep_Hindex
            if (self.DEBUG): print("len(lep_Hindex_vec): ", len(lep_Hindex_vec))
            if len(lep_Hindex_vec)>0:
                for i in range(len(lep_Hindex_vec)):
                    self.branch_values["lep_Hindex"].append(lep_Hindex_vec[i])
            if self.worker.RecoFourMuEvent: finalState = 1
            if self.worker.RecoFourEEvent: finalState = 2
            if self.worker.RecoTwoETwoMuEvent: finalState = 3
            if self.worker.RecoTwoMuTwoEEvent: finalState = 4

        if self.worker.GetZ1_emuCR() and (self.channels == "all"  or self.channels == "2l2v"):
            self.branch_values["foundZZCandidate_2l2nu_emuCR"] = self.worker.ZZSelection_2l2nu()

        if self.DEBUG:
            print("Found ZZ candidate (4l, 2l2q, 2l2nu): ({}, {}, {})".format(self.branch_values["foundZZCandidate_4l"], self.branch_values["foundZZCandidate_2l2q"], self.branch_values["foundZZCandidate_2l2nu"]))

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

        # if (self.branch_definitions["found"]
        if (self.branch_values["foundZZCandidate_4l"] or self.branch_values["foundZZCandidate_2l2q"] or self.branch_values["foundZZCandidate_2l2nu"] or self.branch_values["foundZZCandidate_2l2nu_emuCR"]):
            if self.DEBUG: print("Found ZZ candidate (4l, 2l2q, 2l2nu): ({}, {}, {})".format(self.branch_values["foundZZCandidate_4l"], self.branch_values["foundZZCandidate_2l2q"], self.branch_values["foundZZCandidate_2l2nu"]))
            self.branch_values["pTL1"] = self.worker.pTL1
            self.branch_values["etaL1"] = self.worker.etaL1
            self.branch_values["phiL1"] = self.worker.phiL1
            self.branch_values["massL1"] = self.worker.massL1
            self.branch_values["pTL2"] = self.worker.pTL2
            self.branch_values["etaL2"] = self.worker.etaL2
            self.branch_values["phiL2"] = self.worker.phiL2
            self.branch_values["massL2"] = self.worker.massL2

            if self.branch_values["pTL2"] > self.branch_values["pTL1"]:
                self.branch_values["pTL1"], self.branch_values["pTL2"] = self.branch_values["pTL2"], self.branch_values["pTL1"]
                self.branch_values["etaL1"], self.branch_values["etaL2"] = self.branch_values["etaL2"], self.branch_values["etaL1"]
                self.branch_values["phiL1"], self.branch_values["phiL2"] = self.branch_values["phiL2"], self.branch_values["phiL1"]
                self.branch_values["massL1"], self.branch_values["massL2"] = self.branch_values["massL2"], self.branch_values["massL1"]

            # Kinematics of Z1: Obtained from pair of leptons with mass closest to Z mass
            self.branch_values["pTZ1"] = self.worker.Z1.Pt()
            self.branch_values["etaZ1"] = self.worker.Z1.Eta()
            self.branch_values["phiZ1"] = self.worker.Z1.Phi()
            self.branch_values["massZ1"] = self.worker.Z1.M()

            # Kinematics of Z2: Only for 4l and 2l2q channels
            # For 2l2nu channel, Z2 kinematics are obtained from MET
            # For 2l2q channel, Z2 represents the kinamatics of the boosted Z topology
            self.branch_values["pTZ2"] = self.worker.Z2.Pt()
            self.branch_values["etaZ2"] = self.worker.Z2.Eta()
            self.branch_values["phiZ2"] = self.worker.Z2.Phi()
            self.branch_values["massZ2"] = self.worker.Z2.M()

        if (self.branch_values["foundZZCandidate_2l2q"]):
            keepIt = True
            self.branch_values["foundZZCandidate_2l2q"] = True
            self.passZZ2l2qEvts += 1
            self.CutFlowTable.Fill(4)

            self.branch_values["massZ2_2j"] = self.worker.Z2_2j.M()
            self.branch_values["phiZ2_2j"] = self.worker.Z2_2j.Phi()
            self.branch_values["etaZ2_2j"] = self.worker.Z2_2j.Eta()
            self.branch_values["pTZ2_2j"] = self.worker.Z2_2j.Pt()
            self.branch_values["EneZ2_2j"] = self.worker.Z2_2j.E()

        if (self.branch_values["foundZZCandidate_2l2nu"]):
            keepIt = True
            self.passZZ2l2nuEvts += 1
            self.CutFlowTable.Fill(5)

            #     FatJet_PNZvsQCD = self.worker.FatJet_PNZvsQCD
            #     self.out.fillBranch("FatJet_PNZvsQCD",FatJet_PNZvsQCD)

        if (self.branch_values["foundZZCandidate_2l2nu_emuCR"]):
            keepIt = True
            passZZ2l2nu_emuCR_Selection = True
            self.passZZ2l2nu_emuCR_Evts += 1
            self.CutFlowTable.Fill(6)

        if (self.branch_values["foundZZCandidate_2l2nu"] or self.branch_values["foundZZCandidate_2l2nu_emuCR"]):
            self.branch_values["phiZ2_met"] = self.worker.Z2_met.Phi()
            self.branch_values["pTZ2_met"] = self.worker.Z2_met.Pt()
            self.branch_values["EneZ2_met"] = self.worker.Z2_met.E()
            self.branch_values["MT_2l2nu"] = self.worker.ZZ_metsystem.Mt()

            self.branch_values["HZZ2l2nu_ZZmT"] = self.worker.ZZ_metsystem.Mt()
            self.branch_values["HZZ2l2nu_ZZpT"] = self.worker.ZZ_metsystem.Pt()

            self.branch_values["HZZ2l2nu_ZZmT"] = self.worker.ZZ_metsystem.Mt()
            self.branch_values["HZZ2l2nu_ZZpT"] = self.worker.ZZ_metsystem.Pt()

            #Pz_neutrino = self.worker.Pz_neutrino

            # Define TLorentzVector for VBF jets and get dijet mass
            if HZZ2l2nu_VBFIndexJet1>=0 and HZZ2l2nu_VBFIndexJet2>=0:
                VBF_jet1 = ROOT.TLorentzVector()
                VBF_jet2 = ROOT.TLorentzVector()
                VBF_jet1.SetPtEtaPhiM(jets[HZZ2l2nu_VBFIndexJet1].pt, jets[HZZ2l2nu_VBFIndexJet1].eta, jets[HZZ2l2nu_VBFIndexJet1].phi, jets[HZZ2l2nu_VBFIndexJet1].mass)
                VBF_jet2.SetPtEtaPhiM(jets[HZZ2l2nu_VBFIndexJet2].pt, jets[HZZ2l2nu_VBFIndexJet2].eta, jets[HZZ2l2nu_VBFIndexJet2].phi, jets[HZZ2l2nu_VBFIndexJet2].mass)
                VBF_dijet = VBF_jet1 + VBF_jet2
                if self.DEBUG: print("in .py file: VBF_dijet_mass: ", VBF_dijet.M())

                self.branch_values["HZZ2l2nu_VBFjet1_pT"] = jets[HZZ2l2nu_VBFIndexJet1].pt
                self.branch_values["HZZ2l2nu_VBFjet1_eta"] = jets[HZZ2l2nu_VBFIndexJet1].eta
                self.branch_values["HZZ2l2nu_VBFjet1_phi"] = jets[HZZ2l2nu_VBFIndexJet1].phi
                self.branch_values["HZZ2l2nu_VBFjet1_mass"] = jets[HZZ2l2nu_VBFIndexJet1].mass

                self.branch_values["HZZ2l2nu_VBFjet2_pT"] = jets[HZZ2l2nu_VBFIndexJet2].pt
                self.branch_values["HZZ2l2nu_VBFjet2_eta"] = jets[HZZ2l2nu_VBFIndexJet2].eta
                self.branch_values["HZZ2l2nu_VBFjet2_phi"] = jets[HZZ2l2nu_VBFIndexJet2].phi
                self.branch_values["HZZ2l2nu_VBFjet2_mass"] = jets[HZZ2l2nu_VBFIndexJet2].mass

                self.branch_values["HZZ2l2nu_VBFdijet_mass"] = VBF_dijet.M()
                self.branch_values["HZZ2l2nu_VBFdijet_pT"] = VBF_dijet.Pt()
                self.branch_values["HZZ2l2nu_VBFdijet_E"] = VBF_dijet.E()
                self.branch_values["HZZ2l2nu_VBFdEta_jj"] = abs(VBF_jet1.Eta() - VBF_jet2.Eta())
                self.branch_values["HZZ2l2nu_VBFdPhi_jj"] = abs(VBF_jet1.DeltaPhi(VBF_jet2))
                self.branch_values["HZZ2l2nu_VBFdR_jj"] = VBF_jet1.DeltaR(VBF_jet2)

        if (self.branch_values["foundZZCandidate_4l"]):
            keepIt = True
            self.passZZ4lEvts += 1
            self.CutFlowTable.Fill(3)
            self.branch_values["foundZZCandidate_4l"] = True
            if self.DEBUG: print("Inside 4l loop: ",self.branch_values["foundZZCandidate_4l"])
            self.branch_values["D_CP"] = self.worker.D_CP
            self.branch_values["D_0m"] = self.worker.D_0m
            self.branch_values["D_0hp"] = self.worker.D_0hp
            self.branch_values["D_int"] = self.worker.D_int
            self.branch_values["D_L1"] = self.worker.D_L1
            self.branch_values["D_L1Zg"] = self.worker.D_L1Zg

            self.branch_values["pTL3"] = self.worker.pTL3
            self.branch_values["etaL3"] = self.worker.etaL3
            self.branch_values["phiL3"] = self.worker.phiL3
            self.branch_values["massL3"] = self.worker.massL3
            self.branch_values["pTL4"] = self.worker.pTL4
            self.branch_values["etaL4"] = self.worker.etaL4
            self.branch_values["phiL4"] = self.worker.phiL4
            self.branch_values["massL4"] = self.worker.massL4
            self.branch_values["pTj1"] = self.worker.pTj1
            self.branch_values["etaj1"] = self.worker.etaj1
            self.branch_values["phij1"] = self.worker.phij1
            self.branch_values["massj1"] = self.worker.mj1
            self.branch_values["pTj2"] = self.worker.pTj2
            self.branch_values["etaj2"] = self.worker.etaj2
            self.branch_values["phij2"] = self.worker.phij2
            self.branch_values["massj2"] = self.worker.mj2

            if self.branch_values["pTL4"] > self.branch_values["pTL3"]:
                self.branch_values["pTL3"], self.branch_values["pTL4"] = self.branch_values["pTL4"], self.branch_values["pTL3"]
                self.branch_values["etaL3"], self.branch_values["etaL4"] = self.branch_values["etaL4"], self.branch_values["etaL3"]
                self.branch_values["phiL3"], self.branch_values["phiL4"] = self.branch_values["phiL4"], self.branch_values["phiL3"]
                self.branch_values["massL3"], self.branch_values["massL4"] = self.branch_values["massL4"], self.branch_values["massL3"]

            self.branch_values["pT4l"] = self.worker.ZZsystem.Pt()
            self.branch_values["eta4l"] = self.worker.ZZsystem.Eta()
            self.branch_values["phi4l"] = self.worker.ZZsystem.Phi()
            self.branch_values["mass4l"] = self.worker.ZZsystem.M()
            self.branch_values["rapidity4l"] = self.worker.ZZsystem.Rapidity()
            self.branch_values["njets_pt30_eta4p7"] = self.worker.njets_pt30_eta4p7
            if self.worker.flag4e:
                self.branch_values["mass4e"] = self.branch_values["mass4l"]
            if self.worker.flag2e2mu:
                self.branch_values["mass2e2mu"] = self.branch_values["mass4l"]
            if self.worker.flag4mu:
                self.branch_values["mass4mu"] = self.branch_values["mass4l"]

            if (self.worker.isFSR==False & self.branch_values["passedFullSelection"]):
                self.branch_values["pT4l"] = self.worker.ZZsystemnofsr.Pt()
                self.branch_values["eta4l"] = self.worker.ZZsystemnofsr.Eta()
                self.branch_values["phi4l"] = self.worker.ZZsystemnofsr.Phi()
                self.branch_values["mass4l"] = self.worker.ZZsystemnofsr.M()
                self.branch_values["rapidity4l"] = self.worker.ZZsystemnofsr.Rapidity()

                genWeight = 1 if event.genWeight > 0 else -1
                self.branch_values["Weight"] = genWeight * pileupWeight * dataMCWeight_new * prefiringWeight


        if self.DEBUG:
            print("(found candidates: 2l2q, 2l2nu, 4l): ({:1}, {:1}, {:1}), pTL1: {:7.3f}, pTL2: {:7.3f}, pTZ1: {:7.3f}, pTZ2: {:7.3f}, pTZ2_2j: {:7.3f}, pTZ2_met: {:7.3f}".format(self.branch_values["foundZZCandidate_2l2q"], self.branch_values["foundZZCandidate_2l2nu"], self.branch_values["foundZZCandidate_4l"], self.branch_values["pTL1"], self.branch_values["pTL2"], self.branch_values["pTZ1"], self.branch_values["pTZ2"], self.branch_values["pTZ2_2j"], self.branch_values["pTZ2_met"]))
            if (self.branch_values["foundZZCandidate_2l2q"]):
                print("==> pTL1: {}, \t pTL2: {}".format(self.branch_values["pTL1"], self.branch_values["pTL2"]))

        # Filling the branches with values
        # Fill the branches with the Trigger information for each channel
        for branch, value in self.branch_values.items():
            # Skip filling the branch if it contains "GEN" and the year is not 2022
            # FIXME: Fix it for other years
            if ("GEN" in branch or "lep_genindex" in branch) and self.year != 2022:
                if self.DEBUG: print(f"Skipping branch {branch} for year {self.year}")
                continue

            # if self.DEBUG:
                # print(f"Filling branch {branch} with value {value}")
            if self.DEBUG and "GENlep_MomMomId" in branch:
                print("len(self.branch_values[GENlep_MomMomId]): ", len(self.branch_values["GENlep_MomMomId"]))
            self.out.fillBranch(branch, value)

        return keepIt
