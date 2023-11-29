from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
import ROOT
import yaml
import os
from Helper import *
ROOT.PyConfig.IgnoreCommandLineOptions = True


class HZZAnalysisCppProducer(Module):
    def __init__(self,year,cfgFile,isMC,isFSR):
        base = "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libJHUGenMELAMELA.so" % base)
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libjhugenmela.so" % base)
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libmcfm_707.so" % base)
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libcollier.so" % base)
        if "/H4LTools_cc.so" not in ROOT.gSystem.GetLibraries():
            print("Load C++ module")
            base = "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
            if base:
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/H4LTools.cc+O" % base)
            else:
                base = "$CMSSW_BASE//src/PhysicsTools/NanoAODTools"
                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gROOT.ProcessLine(
                    ".L %s/interface/H4LTools.h" % base)
        if "/RoccoR_cc.so" not in ROOT.gSystem.GetLibraries():
            base = "$CMSSW_BASE//src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
            if base:
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/RoccoR.cc+O" % base)
            else:
                base = "$CMSSW_BASE/src/PhysicsTools/NanoAODTools"
                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gROOT.ProcessLine(
                    ".L %s/interface/RoccoR.h" % base)
        self.year = year
	self.isMC = isMC
        with open(cfgFile, 'r') as ymlfile:
          cfg = yaml.load(ymlfile)
          RoccoRPath = cfg['RoccoRPath']
          self.worker = ROOT.H4LTools(self.year, RoccoRPath)
          self.worker.InitializeElecut(cfg['Electron']['pTcut'],cfg['Electron']['Etacut'],cfg['Electron']['Sip3dcut'],cfg['Electron']['Loosedxycut'],cfg['Electron']['Loosedzcut'],
                                       cfg['Electron']['Isocut'],cfg['Electron']['BDTWP']['LowEta']['LowPT'],cfg['Electron']['BDTWP']['MedEta']['LowPT'],cfg['Electron']['BDTWP']['HighEta']['LowPT'],
                                       cfg['Electron']['BDTWP']['LowEta']['HighPT'],cfg['Electron']['BDTWP']['MedEta']['HighPT'],cfg['Electron']['BDTWP']['HighEta']['HighPT'])
          self.worker.InitializeMucut(cfg['Muon']['pTcut'],cfg['Muon']['Etacut'],cfg['Muon']['Sip3dcut'],cfg['Muon']['Loosedxycut'],cfg['Muon']['Loosedzcut'],cfg['Muon']['Isocut'],
                                       cfg['Muon']['Tightdxycut'],cfg['Muon']['Tightdzcut'],cfg['Muon']['TightTrackerLayercut'],cfg['Muon']['TightpTErrorcut'],cfg['Muon']['HighPtBound'])
          self.worker.InitializeFsrPhotonCut(cfg['FsrPhoton']['pTcut'],cfg['FsrPhoton']['Etacut'],cfg['FsrPhoton']['Isocut'],cfg['FsrPhoton']['dRlcut'],cfg['FsrPhoton']['dRlOverPtcut'])
          self.worker.InitializeJetcut(cfg['Jet']['pTcut'],cfg['Jet']['Etacut'])
          self.worker.InitializeEvtCut(cfg['MZ1cut'],cfg['MZZcut'],cfg['Higgscut']['down'],cfg['Higgscut']['up'],cfg['Zmass'],cfg['MZcut']['down'],cfg['MZcut']['up'])
          self.worker.Initialize2l2qEvtCut(cfg['HZZ2l2q']['Leading_Lep_pT'], cfg['HZZ2l2q']['SubLeading_Lep_pT'], cfg['HZZ2l2q']['Lep_eta'], cfg['HZZ2l2q']['MZLepcut']['down'], cfg['HZZ2l2q']['MZLepcut']['up'])

        self.passtrigEvts = 0
        self.noCutsEvts = 0
        self.passZZ4lEvts = 0
        self.passZZ2l2qEvts = 0
        self.passZZ2l2nuEvts = 0
        MET_sumEt = 0
        self.cfgFile = cfgFile
        self.isMC = isMC
        self.worker.isFSR = isFSR
        pass
    def beginJob(self):
        pass

    def endJob(self):
        print("\n========== Print Cut flow table  ====================\n")
        print("{:27}:{:7} {}".format("Total: ", str(self.noCutsEvts), " Events"))
        print("{:27}:{:7} {}".format("PassTrig: ", str(self.passtrigEvts), " Events"))
        print("{:27}:{:7} {}".format("Pass4eCut: ", str(self.worker.cut4e), " Events"))
        print("Pass4eGhostRemoval: "+str(self.worker.cutghost4e)+" Events")
        print("Pass4eLepPtCut: "+str(self.worker.cutLepPt4e)+" Events")
        print("Pass4eQCDSupress: "+str(self.worker.cutQCD4e)+" Events")
        print("{:27}:{:7} {}".format("PassmZ1mZ2Cut_4e: ", str(self.worker.cutZZ4e), " Events"))
        print("{:27}:{:7} {}".format("Passm4l_105_160_Cut_4e: ", str(self.worker.cutm4l4e), " Events"))
        print("{:27}:{:7} {}".format("Pass4muCut: ", str(self.worker.cut4mu), " Events"))
        print("Pass4muGhostRemoval: "+str(self.worker.cutghost4mu)+" Events")
        print("Pass4muLepPtCut: "+str(self.worker.cutLepPt4mu)+" Events")
        print("Pass4muQCDSupress: "+str(self.worker.cutQCD4mu)+" Events")
        print("{:27}:{:7} {}".format("PassmZ1mZ2Cut_4mu: ", str(self.worker.cutZZ4mu), " Events"))
        print("{:27}:{:7} {}".format("Passm4l_105_160_Cut_4mu: ", str(self.worker.cutm4l4mu), " Events"))
        print("{:27}:{:7} {}".format("Pass2e2muCut: ", str(self.worker.cut2e2mu), " Events"))
        print("Pass2e2muGhostRemoval: "+str(self.worker.cutghost2e2mu)+" Events")
        print("Pass2e2muLepPtCut: "+str(self.worker.cutLepPt2e2mu)+" Events")
        print("Pass2e2muQCDSupress: "+str(self.worker.cutQCD2e2mu)+" Events")
        print("{:27}:{:7} {}".format("PassmZ1mZ2Cut_2e2mu: ", str(self.worker.cutZZ2e2mu), " Events"))
        print("{:27}:{:7} {}".format("Passm4l_105_160_Cut_2e2mu: ", str(self.worker.cutm4l2e2mu), " Events"))
        print("{:27}:{:7} {}".format("PassZZSelection: ", str(self.passZZ4lEvts), " Events"))

        print("\n==================   2l2q    ==============\n")
        print("{:27}:{:7} {}".format("Total: ", str(self.noCutsEvts), " Events"))
        print("{:27}:{:7} {}".format("PassTrig: ", str(self.passtrigEvts), " Events"))
        print("{:27}:{:7} {}".format("Pass2eCut: ", str(self.worker.cut2e), " Events"))
        print("{:27}:{:7} {}".format("Pass2muCut: ", str(self.worker.cut2mu), " Events"))
        print("{:27}:{:7} {}".format("Pass2lCut: ", str(self.worker.cut2l), " Events"))
        print("{:27}:{:7} {}".format("Pass2eCut (40 < mll < 180): ", str(self.worker.cut2e_m40_180), " Events"))
        print("{:27}:{:7} {}".format("Pass2muCut (40 < mll < 180): ", str(self.worker.cut2mu_m40_180), " Events"))
        print("{:27}:{:7} {}".format("Pass2lCut (40 < mll < 180): ", str(self.worker.cut2l_m40_180), " Events"))
        print("{:27}:{:7} {}".format("PassMETcut(gt 150): ", str(self.worker.cutMETlt150), " Events"))
        print("{:27}:{:7} {}".format("Pass2l1JCut: ", str(self.worker.cut2l1J), " Events"))
        print("{:27}:{:7} {}".format("Pass2l2jCut: ", str(self.worker.cut2l2j), " Events"))
        print("{:27}:{:7} {}".format("Pass2l1Jor2jCut: ", str(self.worker.cut2l1Jor2j), " Events"))
        print("{:27}:{:7} {}".format("passZZSelection: ", str(self.passZZ2l2qEvts), " Events"))


        print("\n==================   2l2nu    ==============\n")
        print("{:27}:{:7} {}".format("Total: ", str(self.noCutsEvts), " Events"))
        print("{:27}:{:7} {}".format("PassTrig: ", str(self.passtrigEvts), " Events"))
        print("{:27}:{:7} {}".format("Pass2e_metCut: ", str(self.worker.cut2e_met), " Events"))
        print("{:27}:{:7} {}".format("Pass2mu_metCut: ", str(self.worker.cut2mu_met), " Events"))
        print("{:27}:{:7} {}".format("Pass2l_metCut: ", str(self.worker.cut2l_met), " Events"))
        print("{:27}:{:7} {}".format("Pass2e_metCut (40 < mll < 180): ", str(self.worker.cut2e_met_m40_180), " Events"))
        print("{:27}:{:7} {}".format("Pass2mu_metCut (40 < mll < 180): ", str(self.worker.cut2mu_met_m40_180), " Events"))
        print("{:27}:{:7} {}".format("Pass2l_metCut (40 < mll < 180): ", str(self.worker.cut2l_met_m40_180), " Events"))
        print("{:27}:{:7} {}".format("PassMETcut(gt 150): ", str(self.worker.cutMETgt150), " Events"))
        print("{:27}:{:7} {}".format("Pass2l1metCut: ", str(self.worker.cut2l1met), " Events"))
        print("{:27}:{:7} {}".format("passZZSelection: ", str(self.passZZ2l2nuEvts), " Events"))
        print("\n========== END: Print Cut flow table  ====================\n")

        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree)  # initReaders must be called in beginFile
        self.out = wrappedOutputTree

        # common branches for 4l, 2l2q, 2l2nu channels
        # boolean branch for 4l, 2l2q, 2l2nu channels
        self.out.branch("passZZ4lSelection",  "O")
        self.out.branch("passZZ2l2qSelection",  "O")
        self.out.branch("passZZ2l2nuSelection",  "O")
        self.out.branch("mass4l",  "F")
        self.out.branch("pT4l",  "F")
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
        self.out.branch("massZ2_2j",  "F")
        self.out.branch("phiZ2_2j",  "F")
        self.out.branch("etaZ2_2j",  "F")
        self.out.branch("pTZ2_2j",  "F")
        self.out.branch("EneZ2_2j",  "F")
        self.out.branch("phiZ2_met",  "F")
        self.out.branch("pTZ2_met",  "F")
        self.out.branch("EneZ2_met",  "F")
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
        self.out.branch("mj2",  "F")
        self.out.branch("pTj2",  "F")
        self.out.branch("etaj2",  "F")
        self.out.branch("phij2",  "F")


        self.out.branch("Electron_Fsr_pt",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_eta",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_phi",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Muon_Fsr_pt",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_eta",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_phi",  "F", lenVar = "nMuon_Fsr")

        with open("SyncLepton2018GGH.txt", 'w') as f:
            f.write("Sync data list:"+"\n")

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
        self.noCutsEvts += 1

        massZ2_2j = -999.
        phiZ2_2j = -999.
        etaZ2_2j = -999.
        pTZ2_2j = -999.
        EneZ2_2j = -999.
        phiZ2_met = -999.
        pTZ2_met = -999.
        EneZ2_met = -999.

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

        passedTrig = PassTrig(event, self.cfgFile)
        if (passedTrig==True):
            self.passtrigEvts += 1
        else:
            return keepIt
        electrons = Collection(event, "Electron")
        #electrons = Object(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        jets = Collection(event, "Jet")

        FatJets = Collection(event, "FatJet")
        met = Object(event, "MET", None)
        if isMC:
            genparts = Collection(event, "GenPart")
            for xg in genparts:
                self.worker.SetGenParts(xg.pt)
            for xm in muons:
                self.worker.SetMuonsGen(xm.genPartIdx)

        for xe in electrons:
            self.worker.SetElectrons(xe.pt, xe.eta, xe.phi, xe.mass, xe.dxy,
                                      xe.dz, xe.sip3d, xe.mvaFall17V2Iso, xe.pdgId, xe.pfRelIso03_all)
        for xm in muons:
            self.worker.SetMuons(xm.pt, xm.eta, xm.phi, xm.mass, xm.isGlobal, xm.isTracker,
                                xm.dxy, xm.dz, xm.sip3d, xm.ptErr, xm.nTrackerLayers, xm.isPFcand,
                                 xm.pdgId, xm.charge, xm.pfRelIso03_all)
        for xf in fsrPhotons:
            self.worker.SetFsrPhotons(xf.dROverEt2,xf.eta,xf.phi,xf.pt,xf.relIso03)
        for xj in jets:
            self.worker.SetJets(xj.pt,xj.eta,xj.phi,xj.mass,xj.jetId, xj.btagCSVV2, xj.puId)

        for xj in FatJets:
            self.worker.SetFatJets(xj.pt, xj.eta, xj.phi, xj.msoftdrop, xj.jetId, xj.btagDeepB, xj.particleNet_ZvsQCD)
        self.worker.SetMET(met.pt,met.phi,met.sumEt)
        self.worker.MuonPtCorrection(self.isMC)
        self.worker.LeptonSelection()
        foundZZCandidate_4l = False    # for 4l
        passZZ4lSelection = False
        foundZZCandidate_2l2q = False # for 2l2q
        passZZ2l2qSelection = False
        foundZZCandidate_2l2nu = False # for 2l2nu
        passZZ2l2nuSelection = False
	
	print("=="*51)
        foundZZCandidate_2l2q = self.worker.ZZSelection_2l2q()
        foundZZCandidate_2l2nu = self.worker.ZZSelection_2l2nu()
        foundZZCandidate_4l = self.worker.ZZSelection_4l()

        #if ((self.worker.nTightEle + self.worker.nTightMu == 2) and (not self.worker.nTightMu == 1)):
        #    # This event should belong to either 2l2q or 2l2nu \
        #    # nTightEle + nTightMu == 2 => 2l2q or 2l2nu => (2,0), (0,2), (1,1)
        #    # => Reject (1,1) combination: ( (nTightEle + nTightMu == 2) and (not nTightEle == 1))
        #    # 2nd part is to avoid the situation where we get 1 electron and 1 muon
	#    print("====> met pt: {}, {}, {}".format(met.pt, met.phi, met.sumEt))
        #    foundZZCandidate_2l2q = self.worker.ZZSelection_2l2q()
        #    foundZZCandidate_2l2nu = self.worker.ZZSelection_2l2nu()
        #    print("(2l2q, 2l2nu) = ({}, {})".format(foundZZCandidate_2l2q, foundZZCandidate_2l2nu))
        #    #print("Inside the 2l2q loop: END")

        #elif (self.worker.nTightEle + self.worker.nTightMu >= 4):
        ##if (self.worker.nTightEle + self.worker.nTightMu >= 4):
        #    # This event should belong to 4l; nTightEle + nTightMu >= 4
        #    foundZZCandidate_4l = self.worker.ZZSelection_4l()
        ##if (foundZZCandidate_2l2q and foundZZCandidate_2l2nu ):
        #    #print("both 2l2q and 2l2nu passed the MET selection")
        #    #exit()
        

        if (foundZZCandidate_2l2q):
            keepIt = True
            #print("Hello found 2l2q candidate")
            passZZ2l2qSelection = True
            #print(passZZ2l2qSelection)
            passZZ4lSelection = False
            self.passZZ2l2qEvts += 1
        #     FatJet_PNZvsQCD = self.worker.FatJet_PNZvsQCD
        #     self.out.fillBranch("FatJet_PNZvsQCD",FatJet_PNZvsQCD)
            massZ2_2j = self.worker.Z2_2j.M()  #Anusree
            phiZ2_2j = self.worker.Z2_2j.Phi()
            etaZ2_2j = self.worker.Z2_2j.Eta()
            pTZ2_2j = self.worker.Z2_2j.Pt()
            EneZ2_2j = self.worker.Z2_2j.E()
        
        if (foundZZCandidate_2l2nu):
            keepIt = True
            passZZ2l2nuSelection = True
            self.passZZ2l2nuEvts += 1
        #     FatJet_PNZvsQCD = self.worker.FatJet_PNZvsQCD
        #     self.out.fillBranch("FatJet_PNZvsQCD",FatJet_PNZvsQCD)
            phiZ2_met = self.worker.Z2_met.Phi()
            pTZ2_met = self.worker.Z2_met.Pt()
            EneZ2_met = self.worker.Z2_met.E()
            print("inside 2l2nu loop")
        #self.out.fillBranch("phiZ2_met",phiZ2_met)
        #self.out.fillBranch("pTZ2_met",pTZ2_met)
        #self.out.fillBranch("EneZ2_met",EneZ2_met)
       
        if (foundZZCandidate_4l or foundZZCandidate_2l2q or foundZZCandidate_2l2nu):
            #print("inside loop 4l or 2l2q")
            #print(passZZ2l2qSelection)
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

            pTZ1 = self.worker.Z1.Pt()
            etaZ1 = self.worker.Z1.Eta()
            phiZ1 = self.worker.Z1.Phi()
            massZ1 = self.worker.Z1.M()
            pTZ2 = self.worker.Z2.Pt()
            etaZ2 = self.worker.Z2.Eta()
            phiZ2 = self.worker.Z2.Phi()
            massZ2 = self.worker.Z2.M()

        if (foundZZCandidate_4l):
            keepIt = True
            self.passZZ4lEvts += 1
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
            if self.worker.isFSR==False:
                pT4l = self.worker.ZZsystemnofsr.Pt()
                eta4l = self.worker.ZZsystemnofsr.Eta()
                phi4l = self.worker.ZZsystemnofsr.Phi()
                mass4l = self.worker.ZZsystemnofsr.M()

        self.out.fillBranch("phiZ2_met",phiZ2_met)
        self.out.fillBranch("pTZ2_met",pTZ2_met)
        self.out.fillBranch("EneZ2_met",EneZ2_met)
        
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

        self.out.fillBranch("pTZ1",pTZ1)
        self.out.fillBranch("etaZ1",etaZ1)
        self.out.fillBranch("phiZ1",phiZ1)
        self.out.fillBranch("massZ1",massZ1)
        self.out.fillBranch("pTZ2",pTZ2)
        self.out.fillBranch("etaZ2",etaZ2)
        self.out.fillBranch("phiZ2",phiZ2)
        self.out.fillBranch("massZ2",massZ2)

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

        return keepIt



# define modules using the syntax 'name = lambda : constructor' to avoid
# having them loaded when not needed

#H4LCppModulie() = lambda: HZZAnalysisCppProducer(year)
