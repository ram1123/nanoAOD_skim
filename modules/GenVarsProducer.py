# Reference: https://github.com/latinos/LatinoAnalysis/blob/fae8e13044e23b44961394113a25a8685c4e401b/NanoGardener/python/modules/HiggsGenVarsProducer.py

import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class GenVarsProducer(Module):
    def __init__(self):
        DEBUG = False
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("higgsGenPt", "F")
        self.out.branch("higgsGenEta", "F")
        self.out.branch("higgsGenPhi", "F")
        self.out.branch("higgsGenMass", "F")
        self.out.branch("genV1Pt", "F")
        self.out.branch("genV1Eta", "F")
        self.out.branch("genV1Phi", "F")
        self.out.branch("genV1Mass", "F")
        self.out.branch("genV1DaughterPt", "F", lenVar="nGenV1Daughters")
        self.out.branch("genV1DaughterEta", "F", lenVar="nGenV1Daughters")
        self.out.branch("genV1DaughterPhi", "F", lenVar="nGenV1Daughters")
        self.out.branch("genV1DaughterMass", "F", lenVar="nGenV1Daughters")
        self.out.branch("genV2Pt", "F")
        self.out.branch("genV2Eta", "F")
        self.out.branch("genV2Phi", "F")
        self.out.branch("genV2Mass", "F")
        self.out.branch("genV2DaughterPt", "F", lenVar="nGenV2Daughters")
        self.out.branch("genV2DaughterEta", "F", lenVar="nGenV2Daughters")
        self.out.branch("genV2DaughterPhi", "F", lenVar="nGenV2Daughters")
        self.out.branch("genV2DaughterMass", "F", lenVar="nGenV2Daughters")
        self.out.branch("BoostZ1", "F")
        self.out.branch("BoostZ2", "F")
        self.out.branch("Boostdiff", "F")
        self.out.branch("Pz_neutrino1", "F")
        self.out.branch("Pz_neutrino2", "F")
        #self.out.branch("", "F")
        self.out.branch("delta_pz_neutrino", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def getParentID(self,particle,genParticles):
        if particle.genPartIdxMother == -1: #No parent in record, return ID of original particle
            return particle.pdgId
        elif genParticles[particle.genPartIdxMother].pdgId is particle.pdgId: #'Parent' is self, keep iterating
            return self.getParentID(genParticles[particle.genPartIdxMother],genParticles)
        else: #Found physical parent
            return genParticles[particle.genPartIdxMother].pdgId
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        genParticles = Collection(event, "GenPart")
        ## To obtain GenMET
        genmet = Object(event, "GenMET", None)
        GenMET_pt = None
        GenMET_pt = genmet.pt
        self.DEBUG: print("GenMET_pt: {}".format(GenMET_pt))
        # Loop over gen particles to find Higgs and its each respective decay products. Then keep all kinematics information of Higgs and its respective decay products along with its PDG ID and status flag.

        higgs = None
        v1 = None
        v2 = None
        found_Z1 = False
        found_Z2 = False
        temp_boson = None
        Z1 = ROOT.TLorentzVector(0.0,0.0,0.0, 0.0)
        Z2 = ROOT.TLorentzVector(0.0,0.0,0.0, 0.0)
        v1_decay_products = []
        v2_decay_products = []
        #boost_Z1 = None
        #boost_Z2 = None
        boost_Z1_mag = 0.0
        boost_Z2_mag = 0.0
        boost_diff_mag = 0.0
        neutrino1_pz = 0.0
        neutrino2_pz = 0.0
        delta_pz_neutrino = 0.0
        pz1 = 0.0
        pz2 = 0.0
        E1 = 0.0
        E2 = 0.0
        v1_eta = 0.0

        self.DEBUG: print("length of genParticles: {}".format(len(genParticles)))
        for idx, particle in enumerate(genParticles):
            self.DEBUG: print("DEBUG - line 70: Index: {}, Particle pdgID: {}, Parent ID: {}, MotherIdx: {}, Status: {}".format(idx, particle.pdgId, self.getParentID(particle, genParticles), particle.genPartIdxMother, particle.statusFlags >> 13 & 1))

            if particle.pdgId == 25 and (particle.statusFlags >> 13 & 1):
                higgs = particle
                self.DEBUG: print("DEBUG - line 74 (found higgs): Index: {}, Particle ID: {}, MotherIdx: {}, Parent ID: {}, Status: {}".format(idx, particle.pdgId, self.getParentID(particle, genParticles), particle.genPartIdxMother, particle.statusFlags >> 13 & 1))

            elif (abs(particle.pdgId) == 23) and (particle.statusFlags >> 13 & 1) and self.getParentID(particle, genParticles) == 25:
                self.DEBUG: print("DEBUG - line 76 (found Z boson, daughter of higgs): Index: {}, Particle ID: {}, Parent ID: {}, MotherIdx: {}, Status: {}".format(idx, particle.pdgId, self.getParentID(particle, genParticles), particle.genPartIdxMother, particle.statusFlags >> 13 & 1))

                if v1 is None:
                #if (abs(particle.pdgId) == 23):
                    v1_daughters = []
                    #v2_daughters = []
                    for daughter1 in genParticles:
                        if abs(daughter1.pdgId) in [11, 13] and daughter1.genPartIdxMother == idx and self.getParentID(daughter1, genParticles) == 23 and daughter1.statusFlags >> 13 & 1:
                            v1 = particle
                            #v1_eta = v1.eta
                            #self.out.fillBranch("genV1Eta", v1_eta)
                            found_Z1 = True
                            self.DEBUG: print("DEBUG - line 81 (found Z1 boson, daughter of higgs): Index: {}, Particle ID: {}, Parent ID: {}, MotherIdx: {}, Status: {}".format(idx, v1.pdgId, self.getParentID(v1, genParticles), v1.genPartIdxMother, v1.statusFlags >> 13 & 1))
                            v1_daughters.append(daughter1)
                        #elif temp_boson is None:
                            #temp_boson = particle
                            #if v2 is None:
                                #v2 = temp_boson
                                #v2_daughters.append(daughter1)
                    n = len(v1_daughters)
                    if len(v1_daughters) == 2:
                        for i in range(n):
                             v1_decay_products = v1_daughters
                             self.DEBUG: print("DEBUG - line 92 ( 2 daughters of Z1 boson): Index: {}, Particle ID: {}, Parent ID: {}, MotherIdx: {},  Status: {}".format(idx, v1_decay_products[i].pdgId, self.getParentID(v1_decay_products[i], genParticles), v1_decay_products[i].genPartIdxMother,  v1_decay_products[i].statusFlags >> 13 & 1))

                    #m = len(v2_daughters)
                    #if len(v2_daughters) == 2:
                        #for i in range(m):
                             #v1_decay_products = v1_daughters
                             #self.DEBUG: print("DEBUG - line 108 ( 2 daughters of Z2 boson): Index: {}, Particle ID: {}, Parent ID: {}, MotherIdx: {},  Status: {}".format(idx, v2_decay_products[i].pdgId, self.getParentID(v2_decay_products[i], genParticles), v2_decay_products[i].genPartIdxMother,  v2_decay_products[i].statusFlags >> 13 & 1))
                    #elif abs(v1_daughters[0].pdgId) in [1, 2, 3, 4, 5] and abs(v1_daughters[1].pdgId) in [1, 2, 3, 4, 5]:
                        #v1_decay_products = v1_daughters
                #if v2 is None:
                    #v1 = particle
                    #v2_daughters = []
                    #v1_daughters = []
                if v2 is None:
                    v2_daughters = []
                    for daughter2 in genParticles:
                        if abs(daughter2.pdgId) in [12, 14, 16] and daughter2.genPartIdxMother == idx and self.getParentID(daughter2, genParticles) == 23 and daughter2.statusFlags >> 13 & 1:
                            v2 = particle
                            found_Z2 = True
                            self.DEBUG: print("DEBUG - line 115 (found Z2 boson, daughter of higgs): Index: {}, Particle ID: {}, Parent ID: {}, MotherIdx: {}, Status: {}".format(idx, v2.pdgId, self.getParentID(v2, genParticles), v2.genPartIdxMother, v2.statusFlags >> 13 & 1))
                            v2_daughters.append(daughter2)
                        #elif temp_boson1 is None:
                            #temp_boson1 = particle
                            #if v1 is None:
                                #v1 = temp_boson1
                                #v1_daughters.append(daughter2)
                    m = len(v2_daughters)
                    if len(v2_daughters) == 2:
                        for i in range(m):
                             v2_decay_products = v2_daughters
                             self.DEBUG: print("DEBUG - line 126 ( 2 daughters of Z2 boson): Index: {}, Particle ID: {}, Parent ID: {}, MotherIdx: {},  Status: {}".format(idx, v2_decay_products[i].pdgId, self.getParentID(v2_decay_products[i], genParticles), v2_decay_products[i].genPartIdxMother,  v2_decay_products[i].statusFlags >> 13 & 1))

        if higgs is not None:
            higgs_pt = higgs.pt
            higgs_eta = higgs.eta
            higgs_phi = higgs.phi
            higgs_mass = higgs.mass
        else:
            higgs_pt = -999.
            higgs_eta = -999.
            higgs_phi = -999.
            higgs_mass = -999.

        if v1 is not None:
        #if found_Z1 == True:
            v1_pt = v1.pt
            v1_eta = v1.eta
            v1_phi = v1.phi
            v1_mass = v1.mass

            #Z1 = ROOT.TLorentzVector()
            #Z1.SetPtEtaPhiM(v1_pt, v1_eta, v1_phi, v1_mass)
            #boost_Z1 = Z1.BoostVector()
            #boost_Z1_mag = boost_Z1.Mag()

        else:
            v1_pt = -999.
            v1_eta = -999.
            v1_phi = -999.
            v1_mass = -999.

        if len(v1_decay_products) == 2:
            #v1_decay_products[0] = ROOT.TLorentzVector()
            #v1_decay_products[1] = ROOT.TLorentzVector()
            v1_decay_products_pt = [daughter.pt for daughter in v1_decay_products]
            v1_decay_products_eta = [daughter.eta for daughter in v1_decay_products]
            v1_decay_products_phi = [daughter.phi for daughter in v1_decay_products]
            v1_decay_products_mass = [daughter.mass for daughter in v1_decay_products]
            self.DEBUG: print("v1_decay_products_pt:", v1_decay_products_pt, type(v1_decay_products_pt))
            self.DEBUG: print("v1_decay_products_mass:", v1_decay_products_mass, type(v1_decay_products_mass))
            self.DEBUG: print("v1_decay_products_phi:", v1_decay_products_phi, type(v1_decay_products_phi))
            #pz = v1_decay_products_pt[0] * math.sinh(v1_decay_products_eta[0])
            #self.DEBUG: print("pz of neutrino1:", pz, type(pz))
            #for i in range(2):
                 #v1_decay_products[i] = ROOT.TLorentzVector()
                 #v1_decay_products[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0)
                 #neutrino1_pz = v1_decay_products[0].Pz()
                 #neutrino2_pz = v1_decay_products[1].Pz()
            #self.DEBUG: print("neutrino1_pz:", neutrino1_pz, type(neutrino1_pz))
        else:
            v1_decay_products_pt = [-999.]
            v1_decay_products_eta = [-999.]
            v1_decay_products_phi = [-999.]
            v1_decay_products_mass = [-999.]


        if v2 is not None:
            v2_pt = v2.pt
            v2_eta = v2.eta
            v2_phi = v2.phi
            v2_mass = v2.mass
            self.DEBUG: print("v2_mass:", v2_mass, type(v2_mass))

            #Z2 = ROOT.TLorentzVector()
            #Z2.SetPtEtaPhiM(v2_pt, v2_eta, v2_phi, v2_mass)
            #boost_Z2 = Z2.BoostVector()
            #boost_Z2_mag = boost_Z2.Mag()

        else:
            v2_pt = -999.
            v2_eta = -999.
            v2_phi = -999.
            v2_mass = -999.

        if len(v2_decay_products) == 2:
            v2_decay_products_pt = [daughter.pt for daughter in v2_decay_products]
            v2_decay_products_eta = [daughter.eta for daughter in v2_decay_products]
            v2_decay_products_phi = [daughter.phi for daughter in v2_decay_products]
            v2_decay_products_mass = [daughter.mass for daughter in v2_decay_products]
            self.DEBUG: print("v2_decay_products_pt:", v2_decay_products_pt, type(v2_decay_products_pt))
            pz1 = v2_decay_products_pt[0] * math.sinh(v2_decay_products_eta[0])
            neutrino1_pz = abs(pz1)
            self.DEBUG: print("pz of neutrino1:", neutrino1_pz, type(neutrino1_pz))
            pz2 = v2_decay_products_pt[1] * math.sinh(v2_decay_products_eta[1])
            neutrino2_pz = abs(pz2)
            self.DEBUG: print("pz of neutrino2:", neutrino2_pz, type(neutrino2_pz))
            delta_pz_neutrino = neutrino1_pz - neutrino2_pz
            self.DEBUG: print("delta_pz_neutrino:", delta_pz_neutrino, type(delta_pz_neutrino))

            ###defining Pz of neutrino as a lorentz vector
            v2_decay_products[0] = ROOT.TLorentzVector()
            v2_decay_products[0].SetPtEtaPhiM(v2_decay_products_pt[0], v2_decay_products_eta[0], v2_decay_products_phi[0], v2_decay_products_mass[0])
            v2_decay_products[1] = ROOT.TLorentzVector()
            v2_decay_products[1].SetPtEtaPhiM(v2_decay_products_pt[1], v2_decay_products_eta[1], v2_decay_products_phi[1], v2_decay_products_mass[1])
            #neutrino 1
            pz1 = v2_decay_products[0].Pz()
            self.DEBUG: print("pz of neutrino1_new:", pz1, type(pz1))
            E1 = v2_decay_products[0].E()
            px1 = v2_decay_products[0].Px()
            py1 = v2_decay_products[0].Py()
            p_mag1 = math.sqrt(px1**2 + py1**2 + pz1**2)
            beta1 = v2_decay_products[0].Beta()
            #print("beta 1:", beta1, type(beta1))
            gamma1 = v2_decay_products[0].Gamma()
            #print("gamma1:", gamma1, type(gamma1))
            #pz1_com = gamma1*(pz1 - beta1*E1)
            #print("pz of neutrino1 in COM frame:", pz1_com, type(pz1_com))

            #neutrino 2
            pz2 = v2_decay_products[1].Pz()
            self.DEBUG: print("pz of neutrino2_new:", pz2, type(pz2))
            E2 = v2_decay_products[1].E()
            px2 = v2_decay_products[1].Px()
            py2 = v2_decay_products[1].Py()
        else:
            v2_decay_products_pt = [-999.]
            v2_decay_products_eta = [-999.]
            v2_decay_products_phi = [-999.]
            v2_decay_products_mass = [-999.]
        ##Calculating Pz of neutrino
        #Pz_list = []
        #Pz = ROOT.TMath.Sqrt((v2_mass ** 2) / 4 - GenMET_pt)
        #Pz_list.append(Pz)
        #self.DEBUG: print("Pz:", Pz_list)


        ## Calculating Boost
        if found_Z1 == True and found_Z2 == True:
            Z1.SetPtEtaPhiM(v1_pt, v1_eta, v1_phi, v1_mass)
            boost_Z1 = Z1.BoostVector()
            boost_Z1_mag = boost_Z1.Mag()
            Z2.SetPtEtaPhiM(v2_pt, v2_eta, v2_phi, v2_mass)
            boost_Z2 = Z2.BoostVector()
            boost_Z2_mag = boost_Z2.Mag()
            boost_diff_mag = boost_Z1_mag - boost_Z2_mag
            self.DEBUG: print("delta boost: {}".format(boost_diff_mag))
            #self.out.fillBranch("Boostdiff", boost_diff_mag)

        #self.out.fillBranch("Pz_neutrino", Pz)
        self.out.fillBranch("Pz_neutrino1", pz1)
        self.out.fillBranch("delta_pz_neutrino", delta_pz_neutrino)

        # self.out.fillBranch("Pz_neutrino", Pz)
        self.out.fillBranch("BoostZ1", boost_Z1_mag)
        self.out.fillBranch("BoostZ2", boost_Z2_mag)
        self.out.fillBranch("Boostdiff", boost_diff_mag)
        self.out.fillBranch("higgsGenPt", higgs_pt)
        self.out.fillBranch("higgsGenEta", higgs_eta)
        self.out.fillBranch("higgsGenPhi", higgs_phi)
        self.out.fillBranch("higgsGenMass", higgs_mass)
        self.out.fillBranch("genV1Pt", v1_pt)
        self.out.fillBranch("genV1Eta", v1_eta)
        self.out.fillBranch("genV1Phi", v1_phi)
        self.out.fillBranch("genV1Mass", v1_mass)
        self.out.fillBranch("genV1DaughterPt", v1_decay_products_pt)
        self.out.fillBranch("genV1DaughterEta", v1_decay_products_eta)
        self.out.fillBranch("genV1DaughterPhi", v1_decay_products_phi)
        self.out.fillBranch("genV1DaughterMass", v1_decay_products_mass)
        self.out.fillBranch("genV2Pt", v2_pt)
        self.out.fillBranch("genV2Eta", v2_eta)
        self.out.fillBranch("genV2Phi", v2_phi)
        self.out.fillBranch("genV2Mass", v2_mass)
        self.out.fillBranch("genV2DaughterPt", v2_decay_products_pt)
        self.out.fillBranch("genV2DaughterEta", v2_decay_products_eta)
        self.out.fillBranch("genV2DaughterPhi", v2_decay_products_phi)
        self.out.fillBranch("genV2DaughterMass", v2_decay_products_mass)

        self.DEBUG: print("#######################Event end ################################")
        return True
