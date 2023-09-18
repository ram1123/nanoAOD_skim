from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import yaml
import os
from Helper import *
ROOT.PyConfig.IgnoreCommandLineOptions = True


class HZZAnalysisCppProducer(Module):
    def __init__(self, isMC, year,cfgFile):
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
        self.passtrigEvts = 0
        self.noCutsEvts = 0
        self.passZZEvts = 0
        self.cfgFile = cfgFile
        pass
    def beginJob(self):
        pass

    def endJob(self):
        print("\n========== Print Cut flow table  ====================\n")
        print("{:27}:{:7} {}".format("Total: ", str(self.noCutsEvts), " Events"))
        print("{:27}:{:7} {}".format("PassTrig: ", str(self.passtrigEvts), " Events"))
        print("{:27}:{:7} {}".format("Pass4eCut: ", str(self.worker.cut4e), " Events"))
        print("{:27}:{:7} {}".format("PassmZ1mZ2Cut_4e: ", str(self.worker.cutZZ4e), " Events"))
        print("{:27}:{:7} {}".format("Passm4l_105_160_Cut_4e: ", str(self.worker.cutm4l4e), " Events"))
        print("{:27}:{:7} {}".format("Pass4muCut: ", str(self.worker.cut4mu), " Events"))
        print("{:27}:{:7} {}".format("PassmZ1mZ2Cut_4mu: ", str(self.worker.cutZZ4mu), " Events"))
        print("{:27}:{:7} {}".format("Passm4l_105_160_Cut_4mu: ", str(self.worker.cutm4l4mu), " Events"))
        print("{:27}:{:7} {}".format("Pass2e2muCut: ", str(self.worker.cut2e2mu), " Events"))
        print("{:27}:{:7} {}".format("PassmZ1mZ2Cut_2e2mu: ", str(self.worker.cutZZ2e2mu), " Events"))
        print("{:27}:{:7} {}".format("Passm4l_105_160_Cut_2e2mu: ", str(self.worker.cutm4l2e2mu), " Events"))
        print("{:27}:{:7} {}".format("PassZZSelection: ", str(self.passZZEvts), " Events"))

        print("\n==================   2l2q    ==============\n")
        print("{:27}:{:7} {}".format("Total: ", str(self.noCutsEvts), " Events"))
        print("{:27}:{:7} {}".format("PassTrig: ", str(self.passtrigEvts), " Events"))
        print("{:27}:{:7} {}".format("Pass2eCut: ", str(self.worker.cut2e), " Events"))
        print("{:27}:{:7} {}".format("Pass2muCut: ", str(self.worker.cut2mu), " Events"))
        print("{:27}:{:7} {}".format("Pass2lCut: ", str(self.worker.cut2l), " Events"))
        print("{:27}:{:7} {}".format("Pass2eCut (40 < mll < 180): ", str(self.worker.cut2e_m40_180), " Events"))
        print("{:27}:{:7} {}".format("Pass2muCut (40 < mll < 180): ", str(self.worker.cut2mu_m40_180), " Events"))
        print("{:27}:{:7} {}".format("Pass2lCut (40 < mll < 180): ", str(self.worker.cut2l_m40_180), " Events"))
        print("{:27}:{:7} {}".format("Pass2l1JCut: ", str(self.worker.cut2l1J), " Events"))
        print("{:27}:{:7} {}".format("Pass2l2jCut: ", str(self.worker.cut2l2j), " Events"))
        print("{:27}:{:7} {}".format("Pass2l1Jor2jCut: ", str(self.worker.cut2l1Jor2j), " Events"))
        
        print("\n==================   2l2nu    ==============\n")
        print("{:27}:{:7} {}".format("Total: ", str(self.noCutsEvts), " Events"))
        print("{:27}:{:7} {}".format("PassTrig: ", str(self.passtrigEvts), " Events"))
        print("{:27}:{:7} {}".format("Pass2e_metCut: ", str(self.worker.cut2e_met), " Events"))
        print("{:27}:{:7} {}".format("Pass2mu_metCut: ", str(self.worker.cut2mu_met), " Events"))
        print("{:27}:{:7} {}".format("Pass2l_metCut: ", str(self.worker.cut2l_met), " Events"))
        print("{:27}:{:7} {}".format("Pass2e_metCut (40 < mll < 180): ", str(self.worker.cut2e_met_m40_180), " Events"))
        print("{:27}:{:7} {}".format("Pass2l_metCut (40 < mll < 180): ", str(self.worker.cut2l_met_m40_180), " Events"))
        print("{:27}:{:7} {}".format("Pass2l1metCut: ", str(self.worker.cut2l1met), " Events"))

        print("\n========== END: Print Cut flow table  ====================\n")
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree)  # initReaders must be called in beginFile
        self.out = wrappedOutputTree

        # common branches for 4l, 2l2q, 2l2nu channels
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
        self.out.branch("phiZ2_2l2v",  "F")
        self.out.branch("pTZ2_2l2v",  "F")
        self.out.branch("EneZ2_2l2v",  "F")
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


        self.out.branch("Electron_Fsr_pt",  "F", lenVar = "nElectron")
        self.out.branch("Electron_Fsr_eta",  "F", lenVar = "nElectron")
        self.out.branch("Electron_Fsr_phi",  "F", lenVar = "nElectron")
        self.out.branch("Muon_Fsr_pt",  "F", lenVar = "nMuon")
        self.out.branch("Muon_Fsr_eta",  "F", lenVar = "nMuon")
        self.out.branch("Muon_Fsr_phi",  "F", lenVar = "nMuon")

        # Branches dedicated for 2l2q channel

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
        self.worker.SetObjectNum(event.nElectron,event.nMuon,event.nJet,event.nGenPart,event.nFsrPhoton)

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
        passedTrig = PassTrig(event, self.cfgFile)
        if (passedTrig==True):
            self.passtrigEvts += 1
        else:
            return keepIt
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        jets = Collection(event, "Jet")
        FatJets = Collection(event, "FatJet")
        genparts = Collection(event, "GenPart")
        met = Collection(events, "MET")
        for xe in electrons:
            self.worker.SetElectrons(xe.pt, xe.eta, xe.phi, xe.mass, xe.dxy,
                                      xe.dz, xe.sip3d, xe.mvaFall17V2Iso, xe.pdgId, xe.pfRelIso03_all)
        for xm in muons:
            self.worker.SetMuons(xm.pt, xm.eta, xm.phi, xm.mass, xm.isGlobal, xm.isTracker,
                                xm.dxy, xm.dz, xm.sip3d, xm.ptErr, xm.nTrackerLayers, xm.isPFcand,
                                 xm.pdgId, xm.charge, xm.pfRelIso03_all, xm.genPartIdx)
        for xf in fsrPhotons:
            self.worker.SetFsrPhotons(xf.dROverEt2,xf.eta,xf.phi,xf.pt,xf.relIso03)
        for xj in jets:
            self.worker.SetJets(xj.pt,xj.eta,xj.phi,xj.mass,xj.jetId, xj.btagCSVV2, xj.puId)
        for xj in FatJets:
            self.worker.SetFatJets(xj.pt, xj.eta, xj.phi, xj.msoftdrop, xj.jetId, xj.btagDeepB, xj.particleNet_ZvsQCD)
	for xg in genparts:
            self.worker.SetGenParts(xg.pt)
       # for xn in MET:
        self.worker.SetMET(met.pt,met.phi,met.sumEt)
        self.worker.MuonPtCorrection(self.isMC)
        self.worker.LeptonSelection()
        foundZZCandidate = False    # for 4l
        foundZZCandidate_2l2q = False # for 2l2q
        foundZZCandidate_2l2nu = False # for 2l2nu

        if ((self.worker.nTightEle<2)&(self.worker.nTightMu<2)):
            pass

        if ((self.worker.nTightEle + self.worker.nTightMu == 2) and (not self.worker.nTightMu == 1) and (self.worker.MET_sumEt < 150)):
            # This event should belong to either 2l2q or 2l2nu \
            # nTightEle + nTightMu == 2 => 2l2q or 2l2nu => (2,0), (0,2), (1,1)
            # => Reject (1,1) combination: ( (nTightEle + nTightMu == 2) and (not nTightEle == 1))
            # 2nd part is to avoid the situation where we get 1 electron and 1 muon
            # foundZZCandidate_2l2q = False # print("Inside the 2l2q loop")
            foundZZCandidate_2l2q = self.worker.ZZSelection_2l2q()
            foundZZCandidate_2l2nu = False
            # print("Inside the 2l2q loop: END")
            pass
        elif (self.worker.nTightEle + self.worker.nTightMu >= 4):
            # This event should belong to 4l; nTightEle + nTightMu >= 4
            foundZZCandidate = self.worker.ZZSelection_4l()
        elif ((self.worker.nTightEle + self.worker.nTightMu >= 4) and (self.worker.MET_sumEt > 150)):
            # This event should belong to 4l; nTightEle + nTightMu >= 4
            foundZZCandidate_2l2nu = self.worker.ZZSelection_2l2nu()

        if (foundZZCandidate_2l2q):
            keepIt = True
            self.passZZEvts += 1
        #     FatJet_PNZvsQCD = self.worker.FatJet_PNZvsQCD
        #     self.out.fillBranch("FatJet_PNZvsQCD",FatJet_PNZvsQCD)
            massZ2_2j = self.worker.Z2_2j.M()  
            phiZ2_2j = self.worker.Z2_2j.Phi()
            etaZ2_2j = self.worker.Z2_2j.Eta()
            pTZ2_2j = self.worker.Z2_2j.Pt()
            EneZ2_2j = self.worker.Z2_2j.E()
            self.out.fillBranch("massZ2_2j",massZ2_2j)
            self.out.fillBranch("phiZ2_2j",phiZ2_2j)
            self.out.fillBranch("etaZ2_2j",etaZ2_2j)
            self.out.fillBranch("pTZ2_2j",pTZ2_2j)
            self.out.fillBranch("EneZ2_2j",EneZ2_2j)

        if (foundZZCandidate_2l2nu):
            keepIt = True
            self.passZZEvts += 1
        #     FatJet_PNZvsQCD = self.worker.FatJet_PNZvsQCD
        #     self.out.fillBranch("FatJet_PNZvsQCD",FatJet_PNZvsQCD)
            phiZ2_2l2v = self.worker.Z2_met.Phi()
            pTZ2_2l2v = self.worker.Z2_met.Pt()
            EneZ2_2l2v = self.worker.Z2_met.sumEt()
            self.out.fillBranch("phiZ2_2l2v",phiZ2_2l2v)
            self.out.fillBranch("pTZ2_2l2v",pTZ2_2l2v)
            self.out.fillBranch("EneZ2_2l2v",EneZ2_2l2v)

        if (foundZZCandidate or foundZZCandidate_2l2q or foundZZCandidate_2l2nu):
            keepIt = True
            self.passZZEvts += 1

            pTZ1 = self.worker.Z1.Pt()
            etaZ1 = self.worker.Z1.Eta()
            phiZ1 = self.worker.Z1.Phi()
            massZ1 = self.worker.Z1.M()
            pTZ2 = self.worker.Z2.Pt()
            etaZ2 = self.worker.Z2.Eta()
            phiZ2 = self.worker.Z2.Phi()
            massZ2 = self.worker.Z2.M()

            self.out.fillBranch("pTZ1",pTZ1)
            self.out.fillBranch("etaZ1",etaZ1)
            self.out.fillBranch("phiZ1",phiZ1)
            self.out.fillBranch("massZ1",massZ1)
            self.out.fillBranch("pTZ2",pTZ2)
            self.out.fillBranch("etaZ2",etaZ2)
            self.out.fillBranch("phiZ2",phiZ2)
            self.out.fillBranch("massZ2",massZ2)

        if (foundZZCandidate):
            keepIt = True
            # self.passZZEvts += 1
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

            pT4l = self.worker.ZZsystem.Pt()
            eta4l = self.worker.ZZsystem.Eta()
            phi4l = self.worker.ZZsystem.Phi()
            mass4l = self.worker.ZZsystem.M()
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

        return keepIt


# define modules using the syntax 'name = lambda : constructor' to avoid
# having them loaded when not needed

#H4LCppModule() = lambda: HZZAnalysisCppProducer(year)
