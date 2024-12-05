# Runs Content

| Collection | Description |
| - | - |
| [**LHEPdfSumw**](#lhepdfsumw) | Sum of genEventWeight * LHEPdfWeight[i], divided by genEventSumw |
| [**LHEScaleSumw**](#lhescalesumw) | Sum of genEventWeight * LHEScaleWeight[i], divided by genEventSumw |
| [**genEventCount**](#geneventcount) | event count |
| [**genEventSumw**](#geneventsumw) | sum of gen weights |
| [**genEventSumw2**](#geneventsumw2) | sum of gen (weight^2) |
| [**run**](#run) | run/i |

# Runs detail

## <a id='lhepdfsumw'></a>LHEPdfSumw [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **LHEPdfSumw** | Double_t| Sum of genEventWeight * LHEPdfWeight[i], divided by genEventSumw |
| **nLHEPdfSumw** | UInt_t| Number of entries in LHEPdfSumw |

## <a id='lhescalesumw'></a>LHEScaleSumw [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **LHEScaleSumw** | Double_t| Sum of genEventWeight * LHEScaleWeight[i], divided by genEventSumw |
| **nLHEScaleSumw** | UInt_t| Number of entries in LHEScaleSumw |

## <a id='geneventcount'></a>genEventCount [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **genEventCount** | Long64_t| event count |

## <a id='geneventsumw'></a>genEventSumw [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **genEventSumw** | Double_t| sum of gen weights |

## <a id='geneventsumw2'></a>genEventSumw2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **genEventSumw2** | Double_t| sum of gen (weight^2) |

## <a id='run'></a>run [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **run** | UInt_t| run/i |

# Events Content

| Collection | Description |
| - | - |
| [**ChsMET**](#chsmet) | raw chs PF MET phi |
| [**D**](#d) | MELA discriminant: D_CP |
| [**Electron**](#electron) | Fsr pt of electron |
| [**EneZ2**](#enez2) | Energy of 2nd Z boson obtained from 2 AK4 jets |
| [**HLT**](#hlt) | Trigger/flag bit (process: HLT) |
| [**HZZ2l2nu**](#hzz2l2nu) | VBF topology in 2l2nu channel |
| [**HZZ2l2q**](#hzz2l2q) | PNScore of boosted jet in 2l2q channel |
| [**HZZ2l2qNu**](#hzz2l2qnu) | Electron topology in 2l2q channel |
| [**L1PreFiringWeight**](#l1prefiringweight) | L1 pre-firing event correction weight (1-probability), down var. |
| [**LHEPdfWeight**](#lhepdfweight) | LHE pdf variation weights (w_var / w_nominal) for LHA IDs 306000 - 306102 |
| [**LHEReweightingWeight**](#lhereweightingweight) |  |
| [**LHEScaleWeight**](#lhescaleweight) | LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0  |
| [**LHEWeight**](#lheweight) | Nominal event weight in the LHE file |
| [**MET**](#met) | Delta (METx_mod-METx) Unclustered Energy Up |
| [**MT**](#mt) | Transverse mass of 2l2nu system |
| [**Muon**](#muon) | Fsr pt of muon |
| [**OtherPV**](#otherpv) | Z position of other primary vertices, excluding the main PV |
| [**PSWeight**](#psweight) | PS weights (w_var / w_nominal);   [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5; |
| [**PV**](#pv) | main primary vertex number of degree of freedom |
| [**Pileup**](#pileup) | the true mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled |
| [**PuppiMET**](#puppimet) | phi |
| [**TkMET**](#tkmet) | raw track MET phi |
| [**Weight**](#weight) | Event weight |
| [**btagWeight**](#btagweight) | b-tag event weight for CSVV2 |
| [**eta4l**](#eta4l) | Pseudo-rapidity (eta) of ZZ system |
| [**etaL1**](#etal1) | Pseudo-rapidity (eta) of lepton 1 |
| [**etaL2**](#etal2) | Pseudo-rapidity (eta) of lepton 2 |
| [**etaL3**](#etal3) | Pseudo-rapidity (eta) of lepton 3 |
| [**etaL4**](#etal4) | Pseudo-rapidity (eta) of lepton 4 |
| [**etaZ1**](#etaz1) | Pseudo-rapidity (eta) of Z boson 1 |
| [**etaZ2**](#etaz2) | Pseudo-rapidity (eta) of Z boson 2 |
| [**etaj1**](#etaj1) | Pseudo-rapidity (eta) of leading jet |
| [**etaj2**](#etaj2) | Pseudo-rapidity (eta) of subleading jet |
| [**event**](#event) | event/l |
| [**finalState**](#finalstate) | Final state  identifier for ZZ->4l decay: 1=4mu, 2=4e, 3=2e2mu, 4=2mu2e |
| [**foundZZCandidate**](#foundzzcandidate) | ZZ Candidate found in 4l channel |
| [**genWeight**](#genweight) | generator weight |
| [**isBoosted2l2q**](#isboosted2l2q) | Boosted topology in 2l2q channel |
| [**lep**](#lep) | Reco lepton Higgs indices |
| [**luminosityBlock**](#luminosityblock) | luminosityBlock/i |
| [**mass2e2mu**](#mass2e2mu) | Mass of 2 electrons and 2 muons |
| [**mass4e**](#mass4e) | Mass of 4 electrons |
| [**mass4l**](#mass4l) | Mass of ZZ system in 4l channel |
| [**mass4mu**](#mass4mu) | Mass of 4 muons |
| [**massL1**](#massl1) | Mass of lepton 1 |
| [**massL2**](#massl2) | Mass of lepton 2 |
| [**massL3**](#massl3) | Mass of lepton 3 |
| [**massL4**](#massl4) | Mass of lepton 4 |
| [**massZ1**](#massz1) | Mass of Z boson 1 |
| [**massZ2**](#massz2) | Mass of Z boson 2 |
| [**massj1**](#massj1) | Mass of leading jet |
| [**massj2**](#massj2) | Mass of subleading jet |
| [**njets**](#njets) | Number of jets with pT > 30 GeV and |eta| < 4.7 |
| [**pT4l**](#pt4l) | Transverse momentum (pT) of ZZ system |
| [**pTL1**](#ptl1) | Transverse momentum (pT) of lepton 1 |
| [**pTL2**](#ptl2) | Transverse momentum (pT) of lepton 2 |
| [**pTL3**](#ptl3) | Transverse momentum (pT) of lepton 3 |
| [**pTL4**](#ptl4) | Transverse momentum (pT) of lepton 4 |
| [**pTZ1**](#ptz1) | Transverse momentum (pT) of Z boson 1 |
| [**pTZ2**](#ptz2) | Transverse momentum (pT) of Z boson 2 |
| [**pTj1**](#ptj1) | Transverse momentum (pT) of leading jet |
| [**pTj2**](#ptj2) | Transverse momentum (pT) of subleading jet |
| [**passZZ2l2nu**](#passzz2l2nu) | Pass 2l2nu emu control region selection |
| [**passedFiducialSelection**](#passedfiducialselection) | Passed fiducial selection |
| [**passedFullSelection**](#passedfullselection) | Passed full selection |
| [**phi4l**](#phi4l) | Azimuthal angle (phi) of ZZ system |
| [**phiL1**](#phil1) | Azimuthal angle (phi) of lepton 1 |
| [**phiL2**](#phil2) | Azimuthal angle (phi) of lepton 2 |
| [**phiL3**](#phil3) | Azimuthal angle (phi) of lepton 3 |
| [**phiL4**](#phil4) | Azimuthal angle (phi) of lepton 4 |
| [**phiZ1**](#phiz1) | Azimuthal angle (phi) of Z boson 1 |
| [**phiZ2**](#phiz2) | Azimuthal angle (phi) of Z boson 2 |
| [**phij1**](#phij1) | Azimuthal angle (phi) of leading jet |
| [**phij2**](#phij2) | Azimuthal angle (phi) of subleading jet |
| [**pileupWeight**](#pileupweight) | Pileup reweighting factor |
| [**puWeight**](#puweight) | puWeight/F |
| [**puWeightDown**](#puweightdown) | puWeightDown/F |
| [**puWeightUp**](#puweightup) | puWeightUp/F |
| [**rapidity4l**](#rapidity4l) | Rapidity of ZZ system |
| [**run**](#run) | run/i |

# Events detail

## <a id='chsmet'></a>ChsMET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **ChsMET_phi** | Float_t| raw chs PF MET phi |
| **ChsMET_pt** | Float_t| raw chs PF MET pt |
| **ChsMET_sumEt** | Float_t| raw chs PF scalar sum of Et |

## <a id='d'></a>D [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **D_0hp** | Float_t| MELA discriminant: D_0hp |
| **D_0m** | Float_t| MELA discriminant: D_0m |
| **D_CP** | Float_t| MELA discriminant: D_CP |
| **D_L1** | Float_t| MELA discriminant: D_L1 |
| **D_L1Zg** | Float_t| MELA discriminant: D_L1Zg |
| **D_int** | Float_t| MELA discriminant: D_int |

## <a id='electron'></a>Electron [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Electron_Fsr_eta** | Float_t| Fsr eta of electron |
| **Electron_Fsr_phi** | Float_t| Fsr phi of electron |
| **Electron_Fsr_pt** | Float_t| Fsr pt of electron |
| **nElectron_Fsr** | UInt_t| nElectron_Fsr/i |

## <a id='enez2'></a>EneZ2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **EneZ2_2j** | Float_t| Energy of 2nd Z boson obtained from 2 AK4 jets |
| **EneZ2_met** | Float_t| Energy of 2nd Z boson obtained from MET |

## <a id='hlt'></a>HLT [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **HLT_Dimuon18_PsiPrime_noCorrL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon24_Phi_noCorrL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon24_Upsilon_noCorrL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon25_Jpsi_noCorrL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1** | Bool_t| Trigger/flag bit (process: HLT) |

## <a id='hzz2l2nu'></a>HZZ2l2nu [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **HZZ2l2nu_VBFIndexJet1** | Int_t| Index of leading VBF jet in 2l2nu channel |
| **HZZ2l2nu_VBFIndexJet2** | Int_t| Index of subleading VBF jet in 2l2nu channel |
| **HZZ2l2nu_VBFdEta_jj** | Float_t| Delta Eta between the two VBF jets in 2l2nu channel |
| **HZZ2l2nu_VBFdPhi_jj** | Float_t| Delta Phi between the two VBF jets in 2l2nu channel |
| **HZZ2l2nu_VBFdR_jj** | Float_t| Delta R between the two VBF jets in 2l2nu channel |
| **HZZ2l2nu_VBFdijet_E** | Float_t| Energy of the VBF dijet system in 2l2nu channel |
| **HZZ2l2nu_VBFdijet_mass** | Float_t| Mass of the VBF dijet system in 2l2nu channel |
| **HZZ2l2nu_VBFdijet_pT** | Float_t| Transverse momentum (pT) of the VBF dijet system in 2l2nu channel |
| **HZZ2l2nu_VBFjet1_eta** | Float_t| Pseudo-rapidity (eta) of leading VBF jet in 2l2nu channel |
| **HZZ2l2nu_VBFjet1_mass** | Float_t| Mass of leading VBF jet in 2l2nu channel |
| **HZZ2l2nu_VBFjet1_pT** | Float_t| Transverse momentum (pT) of leading VBF jet in 2l2nu channel |
| **HZZ2l2nu_VBFjet1_phi** | Float_t| Azimuthal angle (phi) of leading VBF jet in 2l2nu channel |
| **HZZ2l2nu_VBFjet2_eta** | Float_t| Pseudo-rapidity (eta) of subleading VBF jet in 2l2nu channel |
| **HZZ2l2nu_VBFjet2_mass** | Float_t| Mass of subleading VBF jet in 2l2nu channel |
| **HZZ2l2nu_VBFjet2_pT** | Float_t| Transverse momentum (pT) of subleading VBF jet in 2l2nu channel |
| **HZZ2l2nu_VBFjet2_phi** | Float_t| Azimuthal angle (phi) of subleading VBF jet in 2l2nu channel |
| **HZZ2l2nu_ZZmT** | Float_t| Transverse mass (mT) of ZZ system in 2l2nu channel |
| **HZZ2l2nu_ZZpT** | Float_t| Transverse momentum (pT) of ZZ system in 2l2nu channel |
| **HZZ2l2nu_ifVBF** | Bool_t| VBF topology in 2l2nu channel |
| **HZZ2l2nu_isEMuCR** | Bool_t| EMU control region in 2l2nu channel |
| **HZZ2l2nu_minDPhi_METAK4** | Float_t| Minimum Delta Phi between MET and AK4 jets in 2l2nu channel |

## <a id='hzz2l2q'></a>HZZ2l2q [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **HZZ2l2q_boostedJet_Index** | Int_t| Index of boosted jet in 2l2q channel |
| **HZZ2l2q_boostedJet_PNScore** | Float_t| PNScore of boosted jet in 2l2q channel |
| **HZZ2l2q_resolvedJet1_Index** | Int_t| Index of 1st resolved jet in 2l2q channel |
| **HZZ2l2q_resolvedJet2_Index** | Int_t| Index of 2nd resolved jet in 2l2q channel |

## <a id='hzz2l2qnu'></a>HZZ2l2qNu [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **HZZ2l2qNu_cutOppositeChargeFlag** | Bool_t| Opposite charge cut passed for 2l2qNu |
| **HZZ2l2qNu_isELE** | Bool_t| Electron topology in 2l2q channel |
| **HZZ2l2qNu_nJets** | Int_t| Number of jets in 2l2qNu channel |
| **HZZ2l2qNu_nLooseBtagJets** | Int_t| Number of loose b-tagged jets in 2l2qNu channel |
| **HZZ2l2qNu_nMediumBtagJets** | Int_t| Number of medium b-tagged jets in 2l2qNu channel |
| **HZZ2l2qNu_nTightBtagJets** | Int_t| Number of tight b-tagged jets in 2l2qNu channel |

## <a id='l1prefiringweight'></a>L1PreFiringWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **L1PreFiringWeight_Dn** | Float_t| L1 pre-firing event correction weight (1-probability), down var. |
| **L1PreFiringWeight_ECAL_Dn** | Float_t| ECAL L1 pre-firing event correction weight (1-probability), down var. |
| **L1PreFiringWeight_ECAL_Nom** | Float_t| ECAL L1 pre-firing event correction weight (1-probability) |
| **L1PreFiringWeight_ECAL_Up** | Float_t| ECAL L1 pre-firing event correction weight (1-probability), up var. |
| **L1PreFiringWeight_Muon_Nom** | Float_t| Muon L1 pre-firing event correction weight (1-probability) |
| **L1PreFiringWeight_Muon_StatDn** | Float_t| Muon L1 pre-firing event correction weight (1-probability), down var. stat. |
| **L1PreFiringWeight_Muon_StatUp** | Float_t| Muon L1 pre-firing event correction weight (1-probability), up var. stat. |
| **L1PreFiringWeight_Muon_SystDn** | Float_t| Muon L1 pre-firing event correction weight (1-probability), down var. syst. |
| **L1PreFiringWeight_Muon_SystUp** | Float_t| Muon L1 pre-firing event correction weight (1-probability), up var. syst. |
| **L1PreFiringWeight_Nom** | Float_t| L1 pre-firing event correction weight (1-probability) |
| **L1PreFiringWeight_Up** | Float_t| L1 pre-firing event correction weight (1-probability), up var. |

## <a id='lhepdfweight'></a>LHEPdfWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **LHEPdfWeight** | Float_t| LHE pdf variation weights (w_var / w_nominal) for LHA IDs 306000 - 306102 |
| **nLHEPdfWeight** | UInt_t|  |

## <a id='lhereweightingweight'></a>LHEReweightingWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **LHEReweightingWeight** | Float_t|  |
| **nLHEReweightingWeight** | UInt_t|  |

## <a id='lhescaleweight'></a>LHEScaleWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **LHEScaleWeight** | Float_t| LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0  |
| **nLHEScaleWeight** | UInt_t|  |

## <a id='lheweight'></a>LHEWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **LHEWeight_originalXWGTUP** | Float_t| Nominal event weight in the LHE file |

## <a id='met'></a>MET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **MET_MetUnclustEnUpDeltaX** | Float_t| Delta (METx_mod-METx) Unclustered Energy Up |
| **MET_MetUnclustEnUpDeltaY** | Float_t| Delta (METy_mod-METy) Unclustered Energy Up |
| **MET_covXX** | Float_t| xx element of met covariance matrix |
| **MET_covXY** | Float_t| xy element of met covariance matrix |
| **MET_covYY** | Float_t| yy element of met covariance matrix |
| **MET_fiducialGenPhi** | Float_t| phi |
| **MET_fiducialGenPt** | Float_t| pt |
| **MET_phi** | Float_t| phi |
| **MET_pt** | Float_t| pt |
| **MET_significance** | Float_t| MET significance |
| **MET_sumEt** | Float_t| scalar sum of Et |
| **MET_sumPtUnclustered** | Float_t| sumPt used for MET significance |

## <a id='mt'></a>MT [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **MT_2l2nu** | Float_t| Transverse mass of 2l2nu system |

## <a id='muon'></a>Muon [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Muon_Fsr_eta** | Float_t| Fsr eta of muon |
| **Muon_Fsr_phi** | Float_t| Fsr phi of muon |
| **Muon_Fsr_pt** | Float_t| Fsr pt of muon |
| **nMuon_Fsr** | UInt_t| nMuon_Fsr/i |

## <a id='otherpv'></a>OtherPV [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **OtherPV_z** | Float_t| Z position of other primary vertices, excluding the main PV |
| **nOtherPV** | UInt_t|  |

## <a id='psweight'></a>PSWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **PSWeight** | Float_t| PS weights (w_var / w_nominal);   [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5; |
| **nPSWeight** | UInt_t|  |

## <a id='pv'></a>PV [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **PV_chi2** | Float_t| main primary vertex reduced chi2 |
| **PV_ndof** | Float_t| main primary vertex number of degree of freedom |
| **PV_npvs** | Int_t| total number of reconstructed primary vertices |
| **PV_npvsGood** | Int_t| number of good reconstructed primary vertices. selection:!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2 |
| **PV_score** | Float_t| main primary vertex score, i.e. sum pt2 of clustered objects |
| **PV_x** | Float_t| main primary vertex position x coordinate |
| **PV_y** | Float_t| main primary vertex position y coordinate |
| **PV_z** | Float_t| main primary vertex position z coordinate |

## <a id='pileup'></a>Pileup [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Pileup_gpudensity** | Float_t| Generator-level PU vertices / mm |
| **Pileup_nPU** | Int_t| the number of pileup interactions that have been added to the event in the current bunch crossing |
| **Pileup_nTrueInt** | Float_t| the true mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled |
| **Pileup_pudensity** | Float_t| PU vertices / mm |
| **Pileup_sumEOOT** | Int_t| number of early out of time pileup |
| **Pileup_sumLOOT** | Int_t| number of late out of time pileup |

## <a id='puppimet'></a>PuppiMET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **PuppiMET_phi** | Float_t| phi |
| **PuppiMET_phiJERDown** | Float_t| JER down phi |
| **PuppiMET_phiJERUp** | Float_t| JER up phi |
| **PuppiMET_phiJESDown** | Float_t| JES down phi |
| **PuppiMET_phiJESUp** | Float_t| JES up phi |
| **PuppiMET_phiUnclusteredDown** | Float_t| Unclustered down phi |
| **PuppiMET_phiUnclusteredUp** | Float_t| Unclustered up phi |
| **PuppiMET_pt** | Float_t| pt |
| **PuppiMET_ptJERDown** | Float_t| JER down pt |
| **PuppiMET_ptJERUp** | Float_t| JER up pt |
| **PuppiMET_ptJESDown** | Float_t| JES down pt |
| **PuppiMET_ptJESUp** | Float_t| JES up pt |
| **PuppiMET_ptUnclusteredDown** | Float_t| Unclustered down pt |
| **PuppiMET_ptUnclusteredUp** | Float_t| Unclustered up pt |
| **PuppiMET_sumEt** | Float_t| scalar sum of Et |

## <a id='tkmet'></a>TkMET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **TkMET_phi** | Float_t| raw track MET phi |
| **TkMET_pt** | Float_t| raw track MET pt |
| **TkMET_sumEt** | Float_t| raw track scalar sum of Et |

## <a id='weight'></a>Weight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Weight** | Float_t| Event weight |

## <a id='btagweight'></a>btagWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **btagWeight_CSVV2** | Float_t| b-tag event weight for CSVV2 |
| **btagWeight_DeepCSVB** | Float_t| b-tag event weight for DeepCSVB |

## <a id='eta4l'></a>eta4l [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **eta4l** | Float_t| Pseudo-rapidity (eta) of ZZ system |

## <a id='etal1'></a>etaL1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **etaL1** | Float_t| Pseudo-rapidity (eta) of lepton 1 |

## <a id='etal2'></a>etaL2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **etaL2** | Float_t| Pseudo-rapidity (eta) of lepton 2 |

## <a id='etal3'></a>etaL3 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **etaL3** | Float_t| Pseudo-rapidity (eta) of lepton 3 |

## <a id='etal4'></a>etaL4 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **etaL4** | Float_t| Pseudo-rapidity (eta) of lepton 4 |

## <a id='etaz1'></a>etaZ1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **etaZ1** | Float_t| Pseudo-rapidity (eta) of Z boson 1 |

## <a id='etaz2'></a>etaZ2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **etaZ2** | Float_t| Pseudo-rapidity (eta) of Z boson 2 |
| **etaZ2_2j** | Float_t| Pseudo-rapidity (eta) of 2nd Z boson obtained from 2 AK4 jets |

## <a id='etaj1'></a>etaj1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **etaj1** | Float_t| Pseudo-rapidity (eta) of leading jet |

## <a id='etaj2'></a>etaj2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **etaj2** | Float_t| Pseudo-rapidity (eta) of subleading jet |

## <a id='event'></a>event [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **event** | ULong64_t| event/l |

## <a id='finalstate'></a>finalState [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **finalState** | Int_t| Final state  identifier for ZZ->4l decay: 1=4mu, 2=4e, 3=2e2mu, 4=2mu2e |

## <a id='foundzzcandidate'></a>foundZZCandidate [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **foundZZCandidate_2l2nu** | Bool_t| ZZ Candidate found in 2l2nu channel |
| **foundZZCandidate_2l2nu_emuCR** | Bool_t| ZZ Candidate found in 2l2nu emu control region |
| **foundZZCandidate_2l2q** | Bool_t| ZZ Candidate found in 2l2q channel |
| **foundZZCandidate_4l** | Bool_t| ZZ Candidate found in 4l channel |

## <a id='genweight'></a>genWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **genWeight** | Float_t| generator weight |

## <a id='isboosted2l2q'></a>isBoosted2l2q [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **isBoosted2l2q** | Bool_t| Boosted topology in 2l2q channel |

## <a id='lep'></a>lep [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **4** | UInt_t| 4/i |
| **lep_Hindex** | Int_t| Reco lepton Higgs indices |

## <a id='luminosityblock'></a>luminosityBlock [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **luminosityBlock** | UInt_t| luminosityBlock/i |

## <a id='mass2e2mu'></a>mass2e2mu [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **mass2e2mu** | Float_t| Mass of 2 electrons and 2 muons |

## <a id='mass4e'></a>mass4e [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **mass4e** | Float_t| Mass of 4 electrons |

## <a id='mass4l'></a>mass4l [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **mass4l** | Float_t| Mass of ZZ system in 4l channel |

## <a id='mass4mu'></a>mass4mu [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **mass4mu** | Float_t| Mass of 4 muons |

## <a id='massl1'></a>massL1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **massL1** | Float_t| Mass of lepton 1 |

## <a id='massl2'></a>massL2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **massL2** | Float_t| Mass of lepton 2 |

## <a id='massl3'></a>massL3 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **massL3** | Float_t| Mass of lepton 3 |

## <a id='massl4'></a>massL4 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **massL4** | Float_t| Mass of lepton 4 |

## <a id='massz1'></a>massZ1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **massZ1** | Float_t| Mass of Z boson 1 |

## <a id='massz2'></a>massZ2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **massZ2** | Float_t| Mass of Z boson 2 |
| **massZ2_2j** | Float_t| Mass of 2nd Z boson obtained from 2 AK4 jets |

## <a id='massj1'></a>massj1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **massj1** | Float_t| Mass of leading jet |

## <a id='massj2'></a>massj2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **massj2** | Float_t| Mass of subleading jet |

## <a id='njets'></a>njets [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **njets_pt30_eta4p7** | Int_t| Number of jets with pT > 30 GeV and |eta| < 4.7 |

## <a id='pt4l'></a>pT4l [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pT4l** | Float_t| Transverse momentum (pT) of ZZ system |

## <a id='ptl1'></a>pTL1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pTL1** | Float_t| Transverse momentum (pT) of lepton 1 |

## <a id='ptl2'></a>pTL2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pTL2** | Float_t| Transverse momentum (pT) of lepton 2 |

## <a id='ptl3'></a>pTL3 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pTL3** | Float_t| Transverse momentum (pT) of lepton 3 |

## <a id='ptl4'></a>pTL4 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pTL4** | Float_t| Transverse momentum (pT) of lepton 4 |

## <a id='ptz1'></a>pTZ1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pTZ1** | Float_t| Transverse momentum (pT) of Z boson 1 |

## <a id='ptz2'></a>pTZ2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pTZ2** | Float_t| Transverse momentum (pT) of Z boson 2 |
| **pTZ2_2j** | Float_t| Transverse momentum (pT) of 2nd Z boson obtained from 2 AK4 jets |
| **pTZ2_met** | Float_t| Transverse momentum (pT) of 2nd Z boson obtained from MET |

## <a id='ptj1'></a>pTj1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pTj1** | Float_t| Transverse momentum (pT) of leading jet |

## <a id='ptj2'></a>pTj2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pTj2** | Float_t| Transverse momentum (pT) of subleading jet |

## <a id='passzz2l2nu'></a>passZZ2l2nu [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **passZZ2l2nu_emuCR_Selection** | Bool_t| Pass 2l2nu emu control region selection |

## <a id='passedfiducialselection'></a>passedFiducialSelection [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **passedFiducialSelection** | Bool_t| Passed fiducial selection |

## <a id='passedfullselection'></a>passedFullSelection [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **passedFullSelection** | Bool_t| Passed full selection |

## <a id='phi4l'></a>phi4l [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **phi4l** | Float_t| Azimuthal angle (phi) of ZZ system |

## <a id='phil1'></a>phiL1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **phiL1** | Float_t| Azimuthal angle (phi) of lepton 1 |

## <a id='phil2'></a>phiL2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **phiL2** | Float_t| Azimuthal angle (phi) of lepton 2 |

## <a id='phil3'></a>phiL3 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **phiL3** | Float_t| Azimuthal angle (phi) of lepton 3 |

## <a id='phil4'></a>phiL4 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **phiL4** | Float_t| Azimuthal angle (phi) of lepton 4 |

## <a id='phiz1'></a>phiZ1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **phiZ1** | Float_t| Azimuthal angle (phi) of Z boson 1 |

## <a id='phiz2'></a>phiZ2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **phiZ2** | Float_t| Azimuthal angle (phi) of Z boson 2 |
| **phiZ2_2j** | Float_t| Azimuthal angle (phi) of 2nd Z boson obtained from 2 AK4 jets |
| **phiZ2_met** | Float_t| Azimuthal angle (phi) of 2nd Z boson obtained from MET |

## <a id='phij1'></a>phij1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **phij1** | Float_t| Azimuthal angle (phi) of leading jet |

## <a id='phij2'></a>phij2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **phij2** | Float_t| Azimuthal angle (phi) of subleading jet |

## <a id='pileupweight'></a>pileupWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **pileupWeight** | Float_t| Pileup reweighting factor |

## <a id='puweight'></a>puWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **puWeight** | Float_t| puWeight/F |

## <a id='puweightdown'></a>puWeightDown [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **puWeightDown** | Float_t| puWeightDown/F |

## <a id='puweightup'></a>puWeightUp [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **puWeightUp** | Float_t| puWeightUp/F |

## <a id='rapidity4l'></a>rapidity4l [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **rapidity4l** | Float_t| Rapidity of ZZ system |

## <a id='run'></a>run [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **run** | UInt_t| run/i |

