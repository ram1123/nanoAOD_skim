from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import ROOT
import yaml
import os
from Helper import *
ROOT.PyConfig.IgnoreCommandLineOptions = True


class HZZAnalysisCppProducer(Module):
    def __init__(self,year,cfgFile,isMC,isFSR,isFiducialAna, DEBUG=False):
        self.loadLibraries()
        self.year = year
        self.isMC = isMC
        self.DEBUG = DEBUG
        self.cfgFile = cfgFile
        self.cfg = self._load_config(cfgFile)
        self.CutFlowTable =  ROOT.TH1F('cutFlow','cutFlow',20, 0, 20)
        self.CutFlowTable.GetXaxis().SetBinLabel(1, "Events")
        self.CutFlowTable.GetXaxis().SetBinLabel(2, "Trigger")
        self.CutFlowTable.GetXaxis().SetBinLabel(3, "4Lepton")
        self.CutFlowTable.GetXaxis().SetBinLabel(4, "4LeptonOSSF")
        self.CutFlowTable.GetXaxis().SetBinLabel(5, "getTightZ")
        self.CutFlowTable.GetXaxis().SetBinLabel(6, "getTightZ1")
        self.CutFlowTable.GetXaxis().SetBinLabel(7, "lep_pTcut")
        self.CutFlowTable.GetXaxis().SetBinLabel(8, "lepdRcut")
        self.CutFlowTable.GetXaxis().SetBinLabel(9, "QCDcut")
        self.CutFlowTable.GetXaxis().SetBinLabel(10, "Smartcut")
        self.CutFlowTable.GetXaxis().SetBinLabel(11, "MZ1MZ2cut")
        self.CutFlowTable.GetXaxis().SetBinLabel(12, "M4Lcut")
        self.CutFlowTable.GetXaxis().SetBinLabel(13, "SR")
        self.CutFlowTable.GetXaxis().SetBinLabel(14, "CR")
        self.CutFlowTable.GetXaxis().SetBinLabel(15, "3Lepton")
        self.CutFlowTable.GetXaxis().SetBinLabel(16, "properID_3lep")
        self.CutFlowTable.GetXaxis().SetBinLabel(17, "3LepDRcut")
        self.CutFlowTable.GetXaxis().SetBinLabel(18, "3LepPtcut")
        self.CutFlowTable.GetXaxis().SetBinLabel(19, "3LepQCDcut")
        self.CutFlowTable.GetXaxis().SetBinLabel(20, "3LepTightZ1cut")

        self.worker = ROOT.H4LTools(self.year, self.isMC)
        self._initialize_worker(self.cfg)
        self.worker.isFSR = isFSR
        self._initialize_counters()
        self.worker.isFiducialAna = isFiducialAna
        pass

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
            fullPath = os.path.join(base_path, 'JHUGenMELA/MELA/data/el9_amd64_gcc12', lib)
            ROOT.gSystem.Load(fullPath)

        # Load the yaml-cpp library
        yaml_cpp_lib_path = os.path.join(yaml_cpp_path, "build")
        ROOT.gSystem.Load(os.path.join(yaml_cpp_lib_path, "libyaml-cpp.so"))

        # Load the C++ module
        if "/H4LTools_cc.so" not in ROOT.gSystem.GetLibraries():
            print("Load H4LTools C++ module")
            if base_path:
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/H4LTools.cc+O" % base_path)
            else:
                base_path = "$CMSSW_BASE//src/PhysicsTools/NanoAODTools"
                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gROOT.ProcessLine(
                    ".L %s/interface/H4LTools.h" % base_path)

    def _load_config(self, cfgFile):
        with open(cfgFile, 'r') as ymlfile:
            return yaml.safe_load(ymlfile)

    def _initialize_worker(self, cfg):
        self.worker.InitializeElecut(cfg['Electron']['pTcut'],cfg['Electron']['Etacut'],cfg['Electron']['Sip3dcut'],cfg['Electron']['Loosedxycut'],cfg['Electron']['Loosedzcut'],
                            cfg['Electron']['Isocut'],cfg['Electron']['BDTWP']['LowEta']['LowPT'],cfg['Electron']['BDTWP']['MedEta']['LowPT'],cfg['Electron']['BDTWP']['HighEta']['LowPT'],
                            cfg['Electron']['BDTWP']['LowEta']['HighPT'],cfg['Electron']['BDTWP']['MedEta']['HighPT'],cfg['Electron']['BDTWP']['HighEta']['HighPT'])
        self.worker.InitializeMucut(cfg['Muon']['pTcut'],cfg['Muon']['Etacut'],cfg['Muon']['Sip3dcut'],cfg['Muon']['Loosedxycut'],cfg['Muon']['Loosedzcut'],cfg['Muon']['Isocut'],
                            cfg['Muon']['Tightdxycut'],cfg['Muon']['Tightdzcut'],cfg['Muon']['TightTrackerLayercut'],cfg['Muon']['TightpTErrorcut'],cfg['Muon']['HighPtBound'])
        self.worker.InitializeFsrPhotonCut(cfg['FsrPhoton']['pTcut'],cfg['FsrPhoton']['Etacut'],cfg['FsrPhoton']['Isocut'],cfg['FsrPhoton']['dRlcut'],cfg['FsrPhoton']['dRlOverPtcut'])
        self.worker.InitializeJetcut(cfg['Jet']['pTcut'],cfg['Jet']['Etacut'])
        self.worker.InitializeEvtCut(cfg['MZ1cut'],cfg['MZZcut'],cfg['Higgscut']['down'],cfg['Higgscut']['up'],cfg['Zmass'],cfg['MZcut']['down'],cfg['MZcut']['up'])

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
        self.passZZEvts = 0

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

        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree)  # initReaders must be called in beginFile
        self.out = wrappedOutputTree
        # Boolean branches for Trigger channels
        for TriggerChannel in self.cfg['TriggerChannels']:
            self.out.branch(TriggerChannel, "O")
        self.out.branch("mass4l",  "F")
        self.out.branch("mass4e",  "F")
        self.out.branch("mass4mu",  "F")
        self.out.branch("mass2e2mu",  "F")
        self.out.branch("pT4l",  "F")
        self.out.branch("rapidity4l",  "F")
        self.out.branch("njets_pt30_eta4p7", "I")
        self.out.branch("nZXCRFailedLeptons", "I")
        self.out.branch("eta4l",  "F")
        self.out.branch("phi4l",  "F")
        self.out.branch("massZ1",  "F")
        self.out.branch("pTZ1",  "F")
        self.out.branch("etaZ1",  "F")
        self.out.branch("phiZ1",  "F")
        self.out.branch("massZ2",  "F")
        self.out.branch("pTZ2",  "F")
        self.out.branch("etaZ2",  "F")
        self.out.branch("phiZ2",  "F")
        self.out.branch("D_CP",  "F")
        self.out.branch("D_0m",  "F")
        self.out.branch("D_0hp",  "F")
        self.out.branch("D_int",  "F")
        self.out.branch("D_L1",  "F")
        self.out.branch("D_L1Zg",  "F")

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
        self.out.branch("mj1",  "F")
        self.out.branch("pTj1",  "F")
        self.out.branch("etaj1",  "F")
        self.out.branch("phij1",  "F")
        self.out.branch("pTj2",  "F")
        self.out.branch("etaj2",  "F")
        self.out.branch("phij2",  "F")
        self.out.branch("mj2",  "F")
        self.out.branch("pileupWeight",  "F")
        self.out.branch("dataMCWeight_new",  "F")
        self.out.branch("prefiringWeight",  "F")
        self.out.branch("passedTrig",  "O")
        self.out.branch("passedFullSelection",  "O")
        self.out.branch("passedZ4lSelection",  "O")
        self.out.branch("passedZ4lZ1LSelection",  "O")
        self.out.branch("passedZ1LSelection",  "O")
        self.out.branch("passedZ4lZXCRSelection",  "O")
        self.out.branch("passedZXCRSelection",  "O")
        self.out.branch("passedFiducialSelection",  "O")
        self.out.branch("CutFlow_4Lepton",  "O")
        self.out.branch("CutFlow_4LeptonOSSF",  "O")
        self.out.branch("CutFlow_getTightZ",  "O")
        self.out.branch("CutFlow_getTightZ1",  "O")
        self.out.branch("CutFlow_lep_pTcut",  "O")
        self.out.branch("CutFlow_lepdRcut",  "O")
        self.out.branch("CutFlow_QCDcut",  "O")
        self.out.branch("CutFlow_Smartcut",  "O")
        self.out.branch("CutFlow_MZ1MZ2cut",  "O")
        self.out.branch("CutFlow_M4Lcut",  "O")
        self.out.branch("CutFlow_SR",  "O")
        self.out.branch("CutFlow_CR",  "O")
        GENHlepNum = 4
        GENZNum = 2
        self.out.branch("lep_Hindex",  "I", lenVar = "GENHlepNum")
        self.out.branch("lep_RelIsoNoFSR",  "F", lenVar = "Lepointer")
        self.out.branch("lep_genindex",  "I", lenVar = "Lepointer")
        self.out.branch("lep_tightId",  "O", lenVar = "Lepointer")
        self.out.branch("lep_id",  "I", lenVar = "Lepointer")
        self.out.branch("lep_pt",  "F", lenVar = "Lepointer")
        self.out.branch("lep_eta",  "F", lenVar = "Lepointer")
        self.out.branch("lep_phi",  "F", lenVar = "Lepointer")
        self.out.branch("lep_mass",  "F", lenVar = "Lepointer")
        self.out.branch("lep_matchedR03_PdgId",  "I", lenVar = "Lepointer")
        self.out.branch("lep_matchedR03_MomId",  "I", lenVar = "Lepointer")
        self.out.branch("lep_matchedR03_MomMomId",  "I", lenVar = "Lepointer")
        self.out.branch("Electron_Fsr_pt",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_eta",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_phi",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Muon_Fsr_pt",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_eta",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_phi",  "F", lenVar = "nMuon_Fsr")

        with open("SyncLepton2018GGH.txt", 'w') as f:
            f.write("Sync data list:"+"\n")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
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
        # do this check at every event, as other modules might have read
        # further branches
        #if event._tree._ttreereaderversion > self._ttreereaderversion:
        #    self.initReaders(event._tree)
        # do NOT access other branches in python between the check/call to
        # initReaders and the call to C++ worker code
        self.worker.Initialize()
        isMC = self.isMC
        self.worker.SetObjectNum(event.nElectron,event.nMuon,event.nJet,event.nFsrPhoton)
        if isMC:
            self.worker.SetObjectNumGen(event.nGenPart)
        keepIt = False

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
        nVECZ = 2
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

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        jets = Collection(event, "Jet")
        if isMC:
            genparts = Collection(event, "GenPart")
            genjets = Collection(event, "GenJet")
            for xg in genparts:
                self.worker.SetGenParts(xg.pt, xg.genPartIdxMother, xg.pdgId)
            for xm in muons:
                self.worker.SetMuonsGen(xm.genPartIdx)
            for xe in electrons:
                self.worker.SetElectronsGen(xe.genPartIdx)
        for xe in electrons:
            self.worker.SetElectrons(xe.pt, xe.eta, xe.phi, xe.mass, xe.dxy,
                                      xe.dz, xe.sip3d, xe.mvaHZZIso, xe.pdgId,xe.charge, xe.pfRelIso03_all)
        for xm in muons:
            self.worker.SetMuons(xm.pt, xm.eta, xm.phi, xm.mass, xm.isGlobal, xm.isTracker,
                                xm.dxy, xm.dz, xm.sip3d, xm.ptErr, xm.nTrackerLayers, xm.nStations, xm.isPFcand,
                                 xm.pdgId, xm.charge, xm.pfRelIso03_all)
        for xf in fsrPhotons:
            self.worker.SetFsrPhotons(xf.dROverEt2,xf.eta,xf.phi,xf.pt,xf.relIso03,xf.electronIdx,xf.muonIdx)
        for xj in jets:
            self.worker.SetJets(xj.pt,xj.eta,xj.phi,xj.mass,xj.jetId, 0.8, 7)
        self.worker.BatchFsrRecovery_Run3()

        self.worker.LeptonSelection()
        foundZCandidate = self.worker.findZCandidate()
        self.worker.findZ1LCandidate()
        if ((self.worker.nTightEle<2)|(self.worker.nTightMu<2)):
            pass

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

        self.worker.findHiggsCandidate()
        foundZZCandidate = self.worker.passedFullSelection
        passedFullSelection = self.worker.passedFullSelection
        passedZ1LSelection = self.worker.passedZ1LSelection
        passedZXCRSelection = self.worker.passedZXCRSelection
        nZXCRFailedLeptons = self.worker.nfailedleptons
        if (passedZ1LSelection): keepIt = True
        if (passedFullSelection): keepIt = True
        if (passedZXCRSelection): keepIt = True
        Lepointer = self.worker.Lepointer
        lep_Hindex = []
        lep_Hindex_vec = self.worker.lep_Hindex
        if len(lep_Hindex_vec)>0:
            for i in range(len(lep_Hindex_vec)):
                lep_Hindex.append(lep_Hindex_vec[i])
        lep_genindex = []
        lep_tightId = []
        lep_tightId_vec = self.worker.lep_tightId
        if len(lep_tightId_vec)>0:
            for i in range(len(lep_tightId_vec)):
                if lep_tightId_vec[i] : lep_tightId.append(1)
                else : lep_tightId.append(0)
        if isMC:
            lep_genindex_vec = self.worker.lep_genindex
            if len(lep_genindex_vec)>0:
                for i in range(len(lep_genindex_vec)):
                    lep_genindex.append(lep_genindex_vec[i])
        lep_pt = []
        lep_eta = []
        lep_phi = []
        lep_mass = []
        lep_id = []
        lep_matchedR03_PdgId = []
        lep_matchedR03_MomId = []
        lep_matchedR03_MomMomId = []
        lep_RelIsoNoFSR = []

        lep_pt_vec = self.worker.lep_pt
        lep_eta_vec = self.worker.lep_eta
        lep_phi_vec = self.worker.lep_phi
        lep_mass_vec = self.worker.lep_mass
        lep_id_vec = self.worker.lep_id
        lep_RelIsoNoFSR_vec = self.worker.lep_RelIsoNoFSR
        lep_matchedR03_PdgId_vec = self.worker.lep_matchedR03_PdgId
        lep_matchedR03_MomId_vec = self.worker.lep_matchedR03_MomId
        lep_matchedR03_MomMomId_vec = self.worker.lep_matchedR03_MomMomId

        if len(lep_pt_vec)>0:
            for i in range(len(lep_pt_vec)):
                lep_pt.append(lep_pt_vec[i])
                lep_eta.append(lep_eta_vec[i])
                lep_phi.append(lep_phi_vec[i])
                lep_mass.append(lep_mass_vec[i])
                lep_id.append(lep_id_vec[i])
                lep_matchedR03_PdgId.append(lep_matchedR03_PdgId_vec[i])
                lep_matchedR03_MomId.append(lep_matchedR03_MomId_vec[i])
                lep_matchedR03_MomMomId.append(lep_matchedR03_MomMomId_vec[i])
                lep_RelIsoNoFSR.append(lep_RelIsoNoFSR_vec[i])

        if (foundZZCandidate):
            self.passZZEvts += 1
        pTZ1 = self.worker.Z1.Pt()
        etaZ1 = self.worker.Z1.Eta()
        phiZ1 = self.worker.Z1.Phi()
        massZ1 = self.worker.Z1.M()
        pTZ2 = self.worker.Z2.Pt()
        etaZ2 = self.worker.Z2.Eta()
        phiZ2 = self.worker.Z2.Phi()
        massZ2 = self.worker.Z2.M()
        D_CP = self.worker.D_CP
        D_0m = self.worker.D_0m
        D_0hp = self.worker.D_0hp
        D_int = self.worker.D_int
        D_L1 = self.worker.D_L1
        D_L1Zg = self.worker.D_L1Zg

        pTL1 = self.worker.pTL1
        etaL1 = self.worker.etaL1
        phiL1 = self.worker.phiL1
        massL1 = self.worker.massL1
        pTL2 = self.worker.pTL2
        etaL2 = self.worker.etaL2
        phiL2 = self.worker.phiL2
        massL2 = self.worker.massL2
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

        CutFlow_4Lepton = self.worker.CutFlow_4Lepton
        CutFlow_4LeptonOSSF = self.worker.CutFlow_4LeptonOSSF
        CutFlow_getTightZ = self.worker.CutFlow_getTightZ
        CutFlow_getTightZ1 = self.worker.CutFlow_getTightZ1
        CutFlow_lep_pTcut = self.worker.CutFlow_lep_pTcut
        CutFlow_lepdRcut = self.worker.CutFlow_lepdRcut
        CutFlow_QCDcut = self.worker.CutFlow_QCDcut
        CutFlow_Smartcut = self.worker.CutFlow_Smartcut
        CutFlow_MZ1MZ2cut = self.worker.CutFlow_MZ1MZ2cut
        CutFlow_M4Lcut = self.worker.CutFlow_M4Lcut
        CutFlow_SR = self.worker.CutFlow_SR
        CutFlow_CR = self.worker.CutFlow_CR
        if CutFlow_4Lepton: self.CutFlowTable.Fill(0)
        if CutFlow_4LeptonOSSF: self.CutFlowTable.Fill(1)
        if CutFlow_getTightZ: self.CutFlowTable.Fill(2)
        if CutFlow_getTightZ1: self.CutFlowTable.Fill(3)
        if CutFlow_lep_pTcut: self.CutFlowTable.Fill(4)
        if CutFlow_lepdRcut: self.CutFlowTable.Fill(5)
        if CutFlow_QCDcut: self.CutFlowTable.Fill(6)
        if CutFlow_Smartcut: self.CutFlowTable.Fill(7)
        if CutFlow_MZ1MZ2cut: self.CutFlowTable.Fill(8)
        if CutFlow_M4Lcut: self.CutFlowTable.Fill(9)
        if CutFlow_SR: self.CutFlowTable.Fill(10)
        if CutFlow_CR: self.CutFlowTable.Fill(11)
        if self.worker.CutFlow_3Lep: self.CutFlowTable.Fill(12)
        if self.worker.CutFlow_properID: self.CutFlowTable.Fill(13)
        if self.worker.CutFlow_3LepDRcut: self.CutFlowTable.Fill(14)
        if self.worker.CutFlow_3LepPtcut: self.CutFlowTable.Fill(15)
        if self.worker.CutFlow_3LepQCDcut: self.CutFlowTable.Fill(16)
        if self.worker.CutFlow_tightZ1cut: self.CutFlowTable.Fill(17)
        if pTL2>pTL1:
            pTL1, pTl2 = pTL2, pTL1
            etaL1, etaL2 = etaL2, etaL1
            phiL1, phiL2 = phiL2, phiL1
            massL1,massL2 = massL2, massL1
        if pTL4>pTL3:
            pTL3, pTL4 = pTL4, pTL3
            etaL3, etaL4 = etaL4, etaL3
            phiL3, phiL4 = phiL4, phiL3
            massL3, massL4 = massL4, massL3
        if (passedFullSelection | passedZXCRSelection):
            pT4l = self.worker.ZZsystem.Pt()
            eta4l = self.worker.ZZsystem.Eta()
            phi4l = self.worker.ZZsystem.Phi()
            mass4l = self.worker.mass4l
            rapidity4l = self.worker.ZZsystem.Rapidity()
        njets_pt30_eta4p7 = self.worker.njets_pt30_eta4p7
        if self.worker.flag4e:
            mass4e = mass4l
        if self.worker.flag2e2mu:
            mass4e = mass4l
        if self.worker.flag4mu:
            mass4mu = mass4l
        if (self.worker.isFSR==False & (passedFullSelection | passedZXCRSelection)):
            pT4l = self.worker.ZZsystemnofsr.Pt()
            eta4l = self.worker.ZZsystemnofsr.Eta()
            phi4l = self.worker.ZZsystemnofsr.Phi()
            mass4l = self.worker.mass4l
            rapidity4l = self.worker.ZZsystemnofsr.Rapidity()


        # Fill the branches with the Trigger information for each channel
        for TriggerChannel in self.cfg['TriggerChannels']:
            self.out.fillBranch(TriggerChannel, TriggerMap[TriggerChannel])
        self.out.fillBranch("mass4l",mass4l)
        self.out.fillBranch("mass4e",mass4e)
        self.out.fillBranch("mass2e2mu",mass2e2mu)
        self.out.fillBranch("mass4mu",mass4mu)
        self.out.fillBranch("pT4l",pT4l)
        self.out.fillBranch("rapidity4l",rapidity4l)
        self.out.fillBranch("njets_pt30_eta4p7",njets_pt30_eta4p7)
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
        self.out.fillBranch("nZXCRFailedLeptons", nZXCRFailedLeptons)
        self.out.fillBranch("passedFullSelection",  passedFullSelection)
        self.out.fillBranch("passedZ4lSelection", passedZ4lSelection)
        self.out.fillBranch("passedZ1LSelection", passedZ1LSelection)
        self.out.fillBranch("passedZ4lZ1LSelection",  passedZ4lZ1LSelection)
        self.out.fillBranch("passedZ4lZXCRSelection",  passedZ4lZXCRSelection)
        self.out.fillBranch("passedZXCRSelection",  passedZXCRSelection)
        self.out.fillBranch("passedFiducialSelection",  passedFiducialSelection)
        self.out.fillBranch("CutFlow_4Lepton",  CutFlow_4Lepton)
        self.out.fillBranch("CutFlow_4LeptonOSSF",  CutFlow_4LeptonOSSF)
        self.out.fillBranch("CutFlow_getTightZ",  CutFlow_getTightZ)
        self.out.fillBranch("CutFlow_getTightZ1",  CutFlow_getTightZ1)
        self.out.fillBranch("CutFlow_lep_pTcut",  CutFlow_lep_pTcut)
        self.out.fillBranch("CutFlow_lepdRcut",  CutFlow_lepdRcut)
        self.out.fillBranch("CutFlow_QCDcut",  CutFlow_QCDcut)
        self.out.fillBranch("CutFlow_Smartcut",  CutFlow_Smartcut)
        self.out.fillBranch("CutFlow_MZ1MZ2cut",  CutFlow_MZ1MZ2cut)
        self.out.fillBranch("CutFlow_M4Lcut",  CutFlow_M4Lcut)
        self.out.fillBranch("CutFlow_SR",  CutFlow_SR)
        self.out.fillBranch("CutFlow_CR",  CutFlow_CR)

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
        self.out.fillBranch("lep_Hindex", lep_Hindex)
        self.out.fillBranch("lep_genindex", lep_genindex)
        self.out.fillBranch("lep_pt", lep_pt)
        self.out.fillBranch("lep_eta", lep_eta)
        self.out.fillBranch("lep_phi", lep_phi)
        self.out.fillBranch("lep_mass", lep_mass)
        self.out.fillBranch("lep_tightId", lep_tightId)
        self.out.fillBranch("lep_id", lep_id)
        self.out.fillBranch("lep_RelIsoNoFSR", lep_RelIsoNoFSR)
        self.out.fillBranch("lep_matchedR03_MomId", lep_matchedR03_MomId)
        self.out.fillBranch("lep_matchedR03_PdgId", lep_matchedR03_PdgId)
        self.out.fillBranch("lep_matchedR03_MomMomId", lep_matchedR03_MomMomId)



        # self.out.fillBranch("nElectron_Fsr", len(electrons))
        # self.out.fillBranch("nMuon_Fsr", len(muons))

        self.out.fillBranch("Electron_Fsr_pt",Electron_Fsr_pt)
        self.out.fillBranch("Electron_Fsr_eta",Electron_Fsr_eta)
        self.out.fillBranch("Electron_Fsr_phi",Electron_Fsr_phi)

        self.out.fillBranch("Muon_Fsr_pt",Muon_Fsr_pt)
        self.out.fillBranch("Muon_Fsr_eta",Muon_Fsr_eta)
        self.out.fillBranch("Muon_Fsr_phi",Muon_Fsr_phi)

        """with open("SyncLepton2018GGH.txt", 'a') as f:
            if(foundZZCandidate):
                f.write(str('%.4f' % event.run)+":"+str('%.4f' % event.luminosityBlock)+":"+str('%.4f' % event.event)+":" \
                        +str('%.4f' % self.worker.pTL1)+":"+str('%.4f' % self.worker.etaL1)+":"+str('%.4f' % self.worker.phiL1)+":"+str('%.4f' % self.worker.massL1)+":" \
                        +str('%.4f' % self.worker.pTL2)+":"+str('%.4f' % self.worker.etaL2)+":"+str('%.4f' % self.worker.phiL2)+":"+str('%.4f' % self.worker.massL2)+":" \
                        +str('%.4f' % self.worker.pTL3)+":"+str('%.4f' % self.worker.etaL3)+":"+str('%.4f' % self.worker.phiL3)+":"+str('%.4f' % self.worker.massL3)+":" \
                        +str('%.4f' % self.worker.pTL4)+":"+str('%.4f' % self.worker.etaL4)+":"+str('%.4f' % self.worker.phiL4)+":"+str('%.4f' % self.worker.massL4)+"\n")
            else:
                f.write(str('%.4f' % event.run)+":"+str('%.4f' % event.luminosityBlock)+":"+str('%.4f' % event.event)+":" \
                        +str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":" \
                        +str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":" \
                        +str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":" \
                        +str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+"\n")"""




        return keepIt


# define modules using the syntax 'name = lambda : constructor' to avoid
# having them loaded when not needed

#H4LCppModule() = lambda: HZZAnalysisCppProducer(year)
