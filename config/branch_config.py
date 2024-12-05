# branch_config.py

# Trigger channel branches
trigger_branches = {
    "TriggerChannels": {
        "type": "O",
        "default": False,
        "title": "Trigger channel: Placeholder for dynamic trigger channels",
    }
}

# Candidate-related branches
candidate_branches = {
    "foundZZCandidate_4l": {"type": "O", "default": False, "title": "ZZ Candidate found in 4l channel"},
    "foundZZCandidate_2l2q": {"type": "O", "default": False, "title": "ZZ Candidate found in 2l2q channel"},
    "foundZZCandidate_2l2nu": {"type": "O", "default": False, "title": "ZZ Candidate found in 2l2nu channel"},
    "foundZZCandidate_2l2nu_emuCR": {"type": "O", "default": False, "title": "ZZ Candidate found in 2l2nu emu control region"},
    "passZZ2l2nu_emuCR_Selection": {"type": "O", "default": False, "title": "Pass 2l2nu emu control region selection"},
    "passedFiducialSelection": {"type": "O", "default": False, "title": "Passed fiducial selection for GEN-level leptons"},
    "isBoosted2l2q": {"type": "O", "default": False, "title": "Boosted topology in 2l2q channel"},
    "HZZ2l2nu_ifVBF": {"type": "O", "default": False, "title": "VBF topology in 2l2nu channel"},
    "HZZ2l2nu_isEMuCR": {"type": "O", "default": False, "title": "EMU control region in 2l2nu channel"},
}

# Lepton-related branches
lepton_branches = {
    f"massL{i}": {"type": "F", "default": -99.0, "title": f"Mass of lepton {i}"} for i in range(1, 5)
} | {
    f"pTL{i}": {"type": "F", "default": -99.0, "title": f"Transverse momentum (pT) of lepton {i}"} for i in range(1, 5)
} | {
    f"etaL{i}": {"type": "F", "default": -99.0, "title": f"Pseudo-rapidity (eta) of lepton {i}"} for i in range(1, 5)
} | {
    f"phiL{i}": {"type": "F", "default": -99.0, "title": f"Azimuthal angle (phi) of lepton {i}"} for i in range(1, 5)
}

# jets branches
jet_branches = {
    "pTj1": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of leading jet"},
    "etaj1": {"type": "F", "default": -99.0, "title": "Pseudo-rapidity (eta) of leading jet"},
    "phij1": {"type": "F", "default": -99.0, "title": "Azimuthal angle (phi) of leading jet"},
    "massj1": {"type": "F", "default": -99.0, "title": "Mass of leading jet"},
    "pTj2": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of subleading jet"},
    "etaj2": {"type": "F", "default": -99.0, "title": "Pseudo-rapidity (eta) of subleading jet"},
    "phij2": {"type": "F", "default": -99.0, "title": "Azimuthal angle (phi) of subleading jet"},
    "massj2": {"type": "F", "default": -99.0, "title": "Mass of subleading jet"}
}

# Z1 and Z2-related branches
z1_branches = {
    "massZ1": {"type": "F", "default": -99.0, "title": "Mass of Z boson 1"},
    "pTZ1": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of Z boson 1"},
    "etaZ1": {"type": "F", "default": -99.0, "title": "Pseudo-rapidity (eta) of Z boson 1"},
    "phiZ1": {"type": "F", "default": -99.0, "title": "Azimuthal angle (phi) of Z boson 1"},
}

# Kinematics of Z2: Only for 4l and 2l2q channels
# For 2l2nu channel, Z2 kinematics are obtained from MET
# For 2l2q channel, Z2 represents the kinamatics of the boosted Z topology only
z2_branches = {
    "massZ2": {"type": "F", "default": -99.0, "title": "Mass of Z boson 2"},
    "pTZ2": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of Z boson 2"},
    "etaZ2": {"type": "F", "default": -99.0, "title": "Pseudo-rapidity (eta) of Z boson 2"},
    "phiZ2": {"type": "F", "default": -99.0, "title": "Azimuthal angle (phi) of Z boson 2"},
}

# ZZ-related branches
zz_branches = {
    "mass4l": {"type": "F", "default": -99.0, "title": "Mass of ZZ system in 4l channel"},
    "GENmass4l": {"type": "F", "default": -99.0, "title": "Generated mass of ZZ system in 4l channel"},
    "mass4e": {"type": "F", "default": -99.0, "title": "Mass of 4 electrons"},
    "mass4mu": {"type": "F", "default": -99.0, "title": "Mass of 4 muons"},
    "mass2e2mu": {"type": "F", "default": -99.0, "title": "Mass of 2 electrons and 2 muons"},
    "pT4l": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of ZZ system"},
    "GENpT4l": {"type": "F", "default": -99.0, "title": "Generated transverse momentum (pT) of ZZ system"},
    "GENrapidity4l": {"type": "F", "default": -99.0, "title": "Generated rapidity of ZZ system"},
    "GENnjets_pt30_eta4p7": {"type": "I", "default": -1, "title": "Number of jets with pT > 30 GeV and |eta| < 4.7"},
    "nGENLeptons": {"type": "I", "default": -1, "title": "Number of GEN leptons"},
    "rapidity4l": {"type": "F", "default": -99.0, "title": "Rapidity of ZZ system"},
    "eta4l": {"type": "F", "default": -99.0, "title": "Pseudo-rapidity (eta) of ZZ system"},
    "phi4l": {"type": "F", "default": -99.0, "title": "Azimuthal angle (phi) of ZZ system"},
    "njets_pt30_eta4p7": {"type": "I", "default": -1, "title": "Number of jets with pT > 30 GeV and |eta| < 4.7"},
    "finalState": {"type": "I", "default": -1, "title": "Final state  identifier for ZZ->4l decay: 1=4mu, 2=4e, 3=2e2mu, 4=2mu2e"},
}

# Channel-specific branches
channel_branches = {
    "HZZ2l2nu_VBFIndexJet1": {"type": "I", "default": -1, "title": "Index of leading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_VBFIndexJet2": {"type": "I", "default": -1, "title": "Index of subleading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_ZZmT": {"type": "F", "default": -99.0, "title": "Transverse mass (mT) of ZZ system in 2l2nu channel"},
    "HZZ2l2nu_ZZpT": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of ZZ system in 2l2nu channel"},
}

TwoL2Q2NuBranches = {
    "HZZ2l2qNu_isELE": {"type": "O", "default": False, "title": "Electron topology in 2l2q channel"},
    "HZZ2l2qNu_cutOppositeChargeFlag": {"type": "O", "default": False, "title": "Opposite charge cut passed for 2l2qNu"},
    "HZZ2l2qNu_nJets": {"type": "I", "default": -1, "title": "Number of jets in 2l2qNu channel"},
    "HZZ2l2qNu_nTightBtagJets": {"type": "I", "default": -1, "title": "Number of tight b-tagged jets in 2l2qNu channel"},
    "HZZ2l2qNu_nMediumBtagJets": {"type": "I", "default": -1, "title": "Number of medium b-tagged jets in 2l2qNu channel"},
    "HZZ2l2qNu_nLooseBtagJets": {"type": "I", "default": -1, "title": "Number of loose b-tagged jets in 2l2qNu channel"},

}

# 2l2nu VBF branches
vbf_2l2nu_branches = {
    "HZZ2l2nu_minDPhi_METAK4": {"type": "F", "default": -99.0, "title": "Minimum Delta Phi between MET and AK4 jets in 2l2nu channel"},
    "HZZ2l2nu_VBFjet1_pT": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of leading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_VBFjet1_eta": {"type": "F", "default": -99.0, "title": "Pseudo-rapidity (eta) of leading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_VBFjet1_phi": {"type": "F", "default": -99.0, "title": "Azimuthal angle (phi) of leading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_VBFjet1_mass": {"type": "F", "default": -99.0, "title": "Mass of leading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_VBFjet2_pT": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of subleading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_VBFjet2_eta": {"type": "F", "default": -99.0, "title": "Pseudo-rapidity (eta) of subleading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_VBFjet2_phi": {"type": "F", "default": -99.0, "title": "Azimuthal angle (phi) of subleading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_VBFjet2_mass": {"type": "F", "default": -99.0, "title": "Mass of subleading VBF jet in 2l2nu channel"},
    "HZZ2l2nu_VBFdijet_mass": {"type": "F", "default": -99.0, "title": "Mass of the VBF dijet system in 2l2nu channel"},
    "HZZ2l2nu_VBFdijet_pT": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of the VBF dijet system in 2l2nu channel"},
    "HZZ2l2nu_VBFdijet_E": {"type": "F", "default": -99.0, "title": "Energy of the VBF dijet system in 2l2nu channel"},
    "HZZ2l2nu_VBFdEta_jj": {"type": "F", "default": -99.0, "title": "Delta Eta between the two VBF jets in 2l2nu channel"},
    "HZZ2l2nu_VBFdPhi_jj": {"type": "F", "default": -99.0, "title": "Delta Phi between the two VBF jets in 2l2nu channel"},
    "HZZ2l2nu_VBFdR_jj": {"type": "F", "default": -99.0, "title": "Delta R between the two VBF jets in 2l2nu channel"},
}

# MELA branches: "D_CP", "D_0m", "D_0hp", "D_int", "D_L1", "D_L1Zg"
mela_branches = {
    "D_CP": {"type": "F", "default": -99.0, "title": "MELA discriminant: D_CP"},
    "D_0m": {"type": "F", "default": -99.0, "title": "MELA discriminant: D_0m"},
    "D_0hp": {"type": "F", "default": -99.0, "title": "MELA discriminant: D_0hp"},
    "D_int": {"type": "F", "default": -99.0, "title": "MELA discriminant: D_int"},
    "D_L1": {"type": "F", "default": -99.0, "title": "MELA discriminant: D_L1"},
    "D_L1Zg": {"type": "F", "default": -99.0, "title": "MELA discriminant: D_L1Zg"},
}

fsr_branches = {
    "Electron_Fsr_pt": {"type": "F", "default": -99.0, "lenVar": "nElectron_Fsr", "title": "Fsr pt of electron"},
    "Electron_Fsr_eta": {"type": "F", "default": -99.0, "lenVar": "nElectron_Fsr", "title": "Fsr eta of electron"},
    "Electron_Fsr_phi": {"type": "F", "default": -99.0, "lenVar": "nElectron_Fsr", "title": "Fsr phi of electron"},
    "Muon_Fsr_pt": {"type": "F", "default": -99.0, "lenVar": "nMuon_Fsr", "title": "Fsr pt of muon"},
    "Muon_Fsr_eta": {"type": "F", "default": -99.0, "lenVar": "nMuon_Fsr", "title": "Fsr eta of muon"},
    "Muon_Fsr_phi": {"type": "F", "default": -99.0, "lenVar": "nMuon_Fsr", "title": "Fsr phi of muon"},
}

# Brancehs for 2l2q channel where Z2 is reconstructed from 1 Fat jet
boostedJetBranches_2l2q = {
    "HZZ2l2q_boostedJet_PNScore": {"type": "F", "default": -99.0, "title": "PNScore of boosted jet in 2l2q channel"},
    "HZZ2l2q_boostedJet_Index": {"type": "I", "default": -1, "title": "Index of boosted jet in 2l2q channel"},
    "HZZ2l2q_resolvedJet1_Index": {"type": "I", "default": -1, "title": "Index of 1st resolved jet in 2l2q channel"},
    "HZZ2l2q_resolvedJet2_Index": {"type": "I", "default": -1, "title": "Index of 2nd resolved jet in 2l2q channel"},
}

# Branches for 2l2q channel where Z2 is reconstructed from 2 jets
TwoLeptonTwoResolvedJetBranches = {
    "massZ2_2j": {"type": "F", "default": -99.0, "title": "Mass of 2nd Z boson obtained from 2 AK4 jets"},
    "phiZ2_2j": {"type": "F", "default": -99.0, "title": "Azimuthal angle (phi) of 2nd Z boson obtained from 2 AK4 jets"},
    "etaZ2_2j": {"type": "F", "default": -99.0, "title": "Pseudo-rapidity (eta) of 2nd Z boson obtained from 2 AK4 jets"},
    "pTZ2_2j": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of 2nd Z boson obtained from 2 AK4 jets"},
    "EneZ2_2j": {"type": "F", "default": -99.0, "title": "Energy of 2nd Z boson obtained from 2 AK4 jets"},
}

# Branches for 2l2v channel where Z2 is reconstructed from MET
TwoLeptonAndMetBranches = {
    "phiZ2_met": {"type": "F", "default": -99.0, "title": "Azimuthal angle (phi) of 2nd Z boson obtained from MET"},
    "pTZ2_met": {"type": "F", "default": -99.0, "title": "Transverse momentum (pT) of 2nd Z boson obtained from MET"},
    "EneZ2_met": {"type": "F", "default": -99.0, "title": "Energy of 2nd Z boson obtained from MET"},
    "MT_2l2nu": {"type": "F", "default": -99.0, "title": "Transverse mass of 2l2nu system"},
}

# Other general branches
general_branches = {
    "EvtNum": {"type": "I", "default": -1, "title": "Event number"},
    "Weight": {"type": "F", "default": 1.0, "title": "Event weight"},
    "pileupWeight": {"type": "F", "default": 1.0, "title": "Pileup reweighting factor"},
    "dataMCWeight_new": {"type": "F", "default": 1.0, "title": "Data/MC reweighting factor"},
}

# Length variables
length_variables = {
    "GENHlepNum": 4,
    "GENZNum": 2
}
# GEN-level branches
gen_level_branches = {
    "GENlep_MomId": {"type": "I", "lenVar": "4", "default": -99, "title": "GEN lepton mother ID"},
    "GENlep_MomMomId": {"type": "I", "lenVar": "4", "default": -99, "title": "GEN lepton grandmother ID"},
    "GENZ_MomId": {"type": "I", "lenVar": "2", "default": -99, "title": "GEN Z boson mother ID"},
    "GENZ_DaughtersId": {"type": "I", "lenVar": "2", "default": -99, "title": "GEN Z boson daughter IDs"},
    "GENlep_Hindex": {"type": "I", "lenVar": "4", "default": -99, "title": "GEN lepton Higgs indices"},
    "lep_Hindex": {"type": "I", "lenVar": "4", "default": -99, "title": "Reco lepton Higgs indices"},
    "GENlep_id": {"type": "I", "lenVar": "4", "default": -99, "title": "GEN lepton IDs"},
    "lep_genindex": {"type": "I", "lenVar": "4", "default": -99, "title": "Reco lepton GEN indices"}
}

# Combine all branches
branch_definitions = {}
branch_definitions.update(trigger_branches)
branch_definitions.update(candidate_branches)
branch_definitions.update(lepton_branches)
branch_definitions.update(jet_branches)
branch_definitions.update(z1_branches)
branch_definitions.update(z2_branches)
branch_definitions.update(zz_branches)
branch_definitions.update(channel_branches)
branch_definitions.update(vbf_2l2nu_branches)
branch_definitions.update(mela_branches)
branch_definitions.update(fsr_branches)
branch_definitions.update(TwoLeptonTwoResolvedJetBranches)
branch_definitions.update(TwoLeptonAndMetBranches)
branch_definitions.update(gen_level_branches)
branch_definitions.update(general_branches)
# branch_definitions.update(length_variables)
branch_definitions.update(TwoL2Q2NuBranches)
branch_definitions.update(boostedJetBranches_2l2q)


if __name__ == "__main__":
    print(branch_definitions)
    print(len(branch_definitions))
