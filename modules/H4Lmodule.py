import ROOT
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import Helper
import Utils

class HZZAnalysisProducer(Module):
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        """self.out.branch("passedTrig",  "I");
        self.out.branch("passedFullSelection",  "I");
        self.out.branch("passedZ4lSelection",  "I");
        self.out.branch("passedQCDcut",  "I");
        self.out.branch("passedZ1LSelection",  "I");
        self.out.branch("passedZ4lZ1LSelection",  "I");
        self.out.branch("passedZ4lZXCRSelection",  "I");
        self.out.branch("passedZXCRSelection",  "I");
        self.out.branch("passedFiducialSelection",  "I");
        self.out.branch("nZXCRFailedLeptons",  "I");"""
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


        """process event, return True (go to next module) or False (fail, go to next event)"""
        #pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """nanoAOD skimming is done considering the final events selection
        for the H4l final state.

        For the H4l semileptonic final state we should have at least 2 tight leptons

        Arguments:
            event {instance of event} -- instance of event

        Returns:
            boolean -- if the event passes skimming then it returns true and
                       go to the next module else returns false and go to
                       the next event.
        """
        elecPtcut = 7
        muPtCut = 5
        tauPtCut = 20
        phoPtCut = 10
        sip3dCut = 4

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
        nZXCRFailedLeptons=0;
        
        HLTs = Collection(event, "HLT")
        if (Helper.PassTrig(HLTs)): passedTrig = True
        else: return keepIt
            
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        taus = Collection(event, "Tau")
        photons = Collection(event, "Photon")

        AllElectrons = Helper.goodLooseElectrons2012(electrons, elecPtcut)
        AllMuons = Helper.goodLooseMuons2012(muons,muPtCut)
        AllTaus =  Helper.goodLooseTaus2012(taus,tauPtCut)
        AllPhotons = Helper.goodLoosePhotons2015(photons,phoPtCut)
        
        #loose leptons
        recoMuons = Helper.goodMuons2015_noIso_noPf(AllMuons,muPtCut,sip3dCut)
        recoElectrons = Helper.goodElectrons2015_noIso_noBdt(AllElectrons, elecPtcut, sip3dCut)
        recoTaus = Helper.goodTaus2015(AllTaus, tauPtCut)
        recoPhotons = Helper.goodPhotons2015(AllPhotons, phoPtCut)

        Ele_tight_ID = Helper.passTight_BDT_Id(recoElectrons)
        Mu_tight_ID = Helper.passTight_Id(recoMuons)

        Z1 = ROOT.TLorentzVector()
        Z2 = ROOT.TLorentzVector()

        foundZZCandidate, Z1, Z2 = Utils.ZZSelection(recoElectrons, recoMuons, Ele_tight_ID, Mu_tight_ID)
        
        if (foundZZCandidate):
            keepIt = True
            pTZ1 = Z1.Pt()
            etaZ1 = Z1.Eta()
            phiZ1 = Z1.Phi()
            massZ1 = Z1.M()
            pTZ2 = Z2.Pt()
            etaZ2 = Z2.Eta()
            phiZ2 = Z2.Phi()
            massZ2 = Z2.M()

            ZZsystem = ROOT.TLorentzVector()
            ZZsystem = Z1 + Z2
            pT4l = ZZsystem.Pt()
            eta4l = ZZsystem.Eta()
            phi4l = ZZsystem.Phi()
            mass4l = ZZsystem.M()
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
        """self.out.fillBranch("passedTrig",passedTrig)
        self.out.fillBranch("passedFullSelection",passedFullSelection)
        self.out.fillBranch("passedZ4lSelection",passedZ4lSelection)
        self.out.fillBranch("passedQCDcut",passedQCDcut)
        self.out.fillBranch("passedZ1LSelection",passedZ1LSelection)
        self.out.fillBranch("passedZ4lZ1LSelection",passedZ4lZ1LSelection)
        self.out.fillBranch("passedZ4lZXCRSelection",passedZ4lZXCRSelection)
        self.out.fillBranch("passedZXCRSelection",passedZXCRSelection)
        self.out.fillBranch("passedFiducialSelection",passedFiducialSelection)
        self.out.fillBranch("nZXCRFailedLeptons",nZXCRFailedLeptons)"""

        return keepIt


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
HZZAnalysisModule = lambda : HZZAnalysisProducer() #(jetSelection= lambda j : j.pt > 30)