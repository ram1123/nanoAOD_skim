from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import yaml
import os
from Helper import *
ROOT.PyConfig.IgnoreCommandLineOptions = True


class HZZAnalysisCppProducer(Module):
    def __init__(self,year,cfgFile,isMC,isFSR,isFiducialAna):
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
        self.year = year
        with open(cfgFile, 'r') as ymlfile:
          cfg = yaml.load(ymlfile)
          self.worker = ROOT.H4LTools(self.year)
          self.worker.InitializeElecut(cfg['Electron']['pTcut'],cfg['Electron']['Etacut'],cfg['Electron']['Sip3dcut'],cfg['Electron']['Loosedxycut'],cfg['Electron']['Loosedzcut'],
                                       cfg['Electron']['Isocut'],cfg['Electron']['BDTWP']['LowEta']['LowPT'],cfg['Electron']['BDTWP']['MedEta']['LowPT'],cfg['Electron']['BDTWP']['HighEta']['LowPT'],
                                       cfg['Electron']['BDTWP']['LowEta']['HighPT'],cfg['Electron']['BDTWP']['MedEta']['HighPT'],cfg['Electron']['BDTWP']['HighEta']['HighPT'])
          self.worker.InitializeMucut(cfg['Muon']['pTcut'],cfg['Muon']['Etacut'],cfg['Muon']['Sip3dcut'],cfg['Muon']['Loosedxycut'],cfg['Muon']['Loosedzcut'],cfg['Muon']['Isocut'],
                                       cfg['Muon']['Tightdxycut'],cfg['Muon']['Tightdzcut'],cfg['Muon']['TightTrackerLayercut'],cfg['Muon']['TightpTErrorcut'],cfg['Muon']['HighPtBound'])
          self.worker.InitializeFsrPhotonCut(cfg['FsrPhoton']['pTcut'],cfg['FsrPhoton']['Etacut'],cfg['FsrPhoton']['Isocut'],cfg['FsrPhoton']['dRlcut'],cfg['FsrPhoton']['dRlOverPtcut'])
          self.worker.InitializeJetcut(cfg['Jet']['pTcut'],cfg['Jet']['Etacut'])
          self.worker.InitializeEvtCut(cfg['MZ1cut'],cfg['MZZcut'],cfg['Higgscut']['down'],cfg['Higgscut']['up'],cfg['Zmass'],cfg['MZcut']['down'],cfg['MZcut']['up'])

        self.passtrigEvts = 0
        self.passZZEvts = 0
        self.cfgFile = cfgFile
        self.isMC = isMC
        self.worker.isFSR = isFSR
        self.worker.isFiducialAna = isFiducialAna
        pass
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
        passedTrig = PassTrig(event, self.cfgFile)
        if (passedTrig==True):
            self.passtrigEvts += 1
        else:
            return keepIt
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        jets = Collection(event, "Jet")
        if isMC:
            genparts = Collection(event, "GenPart")
            for xg in genparts:
                self.worker.SetGenParts(xg.pt)
            for xm in muons:
                self.worker.SetMuonsGen(xm.genPartIdx)
        for xe in electrons:
            self.worker.SetElectrons(xe.pt, xe.eta, xe.phi, xe.mass, xe.dxy,
                                      xe.dz, xe.sip3d, xe.mvaIso_Fall17V2, xe.pdgId, xe.pfRelIso03_all)
        for xm in muons:
            self.worker.SetMuons(xm.corrected_pt, xm.eta, xm.phi, xm.mass, xm.isGlobal, xm.isTracker,
                                xm.dxy, xm.dz, xm.sip3d, xm.ptErr, xm.nTrackerLayers, xm.isPFcand,
                                 xm.pdgId, xm.charge, xm.pfRelIso03_all)
        for xf in fsrPhotons:
            self.worker.SetFsrPhotons(xf.dROverEt2,xf.eta,xf.phi,xf.pt,xf.relIso03)
        for xj in jets:
            self.worker.SetJets(xj.pt,xj.eta,xj.phi,xj.mass,xj.jetId, xj.btagCSVV2, 7)

        self.worker.LeptonSelection()
        foundZCandidate = self.worker.findZCandidate()
        self.worker.findZ1LCandidate()
        if ((self.worker.nTightEle<2)&(self.worker.nTightMu<2)):
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

        foundZZCandidate = self.worker.ZZSelection()

        if (foundZZCandidate):
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


            pT4l = self.worker.ZZsystem.Pt()
            eta4l = self.worker.ZZsystem.Eta()
            phi4l = self.worker.ZZsystem.Phi()
            mass4l = self.worker.ZZsystem.M()
            if self.worker.isFSR==False:
                pT4l = self.worker.ZZsystemnofsr.Pt()
                eta4l = self.worker.ZZsystemnofsr.Eta()
                phi4l = self.worker.ZZsystemnofsr.Phi()
                mass4l = self.worker.ZZsystemnofsr.M()
            self.out.fillBranch("mass4l",mass4l)
            self.out.fillBranch("pT4l",pT4l)
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
