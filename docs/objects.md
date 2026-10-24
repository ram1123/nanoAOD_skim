# Object selection

Thresholds come from `config/Input_<year>.yml`, which the Python module pushes
into the C++ worker through `Initialize*` setters. The logic is in
`src/H4LTools.cc`. Object cuts are **identical in all four years**; only the
triggers and the b-tag working points change (see
[Event selection](event-selection.md#triggers)).

## Muons

Muon p<sub>T</sub> is the **Rochester-corrected** `Muon_corrected_pt`, from
`muonScaleRes20XX`, which runs before the analysis module.

| step | requirement | YAML key |
| --- | --- | --- |
| kinematics | p<sub>T</sub> > 25 GeV, \|&eta;\| < 2.4 | `Muon.pTcut`, `Muon.Etacut` |
| ID | `Muon_mediumId` and (`isGlobal` or `isTracker`) | |
| impact parameter | \|d<sub>xy</sub>\| < 0.045 cm, \|d<sub>z</sub>\| < 0.2 cm | `Loosedxycut`, `Loosedzcut` |
| "tight" flag | p<sub>T</sub> < 200: `isPFcand`. p<sub>T</sub> &ge; 200: `isPFcand` **or** (&sigma;<sub>p<sub>T</sub></sub>/p<sub>T</sub> < 0.3, \|d<sub>xy</sub>\| < 0.045, \|d<sub>z</sub>\| < 0.2, tracker layers > 5) | `HighPtBound`, `TightpTErrorcut`, `TightTrackerLayercut` |
| isolation | `pfRelIso03_all` &minus; p<sub>T</sub><sup>FSR&gamma;</sup>/p<sub>T</sub> < 0.20 | `Muon.Isocut` |

A muon that passes all rows is a **tight muon**. The cut flow's `cut_mu_*` rows
count events with at least 2 muons passing these cuts cumulatively.

## Electrons

| step | requirement | YAML key |
| --- | --- | --- |
| kinematics | p<sub>T</sub> > 25 GeV, \|&eta;\| < 2.5 | `Electron.pTcut`, `Electron.Etacut` |
| ECAL gap | reject 1.4442 < \|&eta;<sub>SC</sub>\| < 1.5660, with &eta;<sub>SC</sub> = &eta; + `deltaEtaSC` | |
| impact parameter | \|d<sub>xy</sub>\| < 0.045 cm, \|d<sub>z</sub>\| < 0.2 cm | `Loosedxycut`, `Loosedzcut` |
| ID | `Electron_mvaIso_WP90` (v15; v9 used `mvaFall17V2Iso_WP90`) | |
| isolation | `pfRelIso03_all` &minus; p<sub>T</sub><sup>FSR&gamma;</sup>/p<sub>T</sub> < 0.15 | `Electron.Isocut` |

Electron&ndash;muon cross-cleaning (`FlagElectronMuonOverlap`, &Delta;R < 0.05 to a
loose muon) is implemented but **disabled**, because AN-2016/325 does not use it.
The `MuonCleaning` YAML block is read but has no effect.

## FSR photons

FSR recovery is on (`isFSR = True`). It uses `FsrPhoton` candidates with
p<sub>T</sub> > 2 GeV, \|&eta;\| < 2.4 and `relIso03` < 0.8. A candidate must also
satisfy &Delta;R(&ell;, &gamma;) < 0.5 and `dROverEt2` < 0.012. For each lepton,
the candidate with the smallest `dROverEt2` is added to the lepton four-vector and
subtracted from the lepton's isolation. Z candidates are built from the
FSR-corrected leptons. The `Electron_Fsr_*` and `Muon_Fsr_*` branches store the
corrected kinematics.

## Veto leptons (extra-lepton veto)

Leptons beyond the Z&#8321; or e&mu; pair count into `HZZ2l2nu_nExtraLep`. These
are AN-2016/325 &sect;4.4.4 values, not re-verified for UL / v15.

| object | definition | YAML (`HZZ2l2nu.VetoLep`) |
| --- | --- | --- |
| veto electron | p<sub>T</sub> > 10, \|&eta;\| < 2.5, `cutBased` &ge; 2 (Fall17V2 Loose), `pfRelIso03_all` < 0.15 | `Ele_pT`, `Ele_cutBased`, `Ele_Iso` |
| veto muon | p<sub>T</sub> > 10, \|&eta;\| < 2.4, `isPFcand` and (`isGlobal` or `isTracker`), `pfRelIso03_all` < 0.2 | `Mu_pT`, `Mu_Iso` |
| any further tight lepton | counts even if it fails the veto ID, so the event has **exactly two** tight leptons | |
| soft muon | `Muon_softId`, p<sub>T</sub> > 3, \|&eta;\| < 2.4, not already a veto muon. Two counts: `HZZ2l2nu_nExtraSoftMu` (AN definition) and `HZZ2l2nu_nExtraSoftMuPV` (adds \|d<sub>z</sub>\| < 0.5, \|d<sub>xy</sub>\| < 0.2) | `SoftMu_pT` (< 0: off) |

## Photons (&gamma;+jets CR only, `--channels 1g`) { #photons }

| requirement | YAML (`GammaCR`) |
| --- | --- |
| p<sub>T</sub> > 55 GeV | `pTcut` |
| barrel supercluster (`Photon_isScEtaEB`) | |
| `Photon_cutBased` &ge; 3 (tight) | `IDcut` |
| `electronVeto` and no `pixelSeed` | |
| **exactly one** such photon in the event | |

## Jets

v15 jets are **AK4 PUPPI** and carry no stored jet ID or PU ID.

| step | requirement | YAML key |
| --- | --- | --- |
| kinematics | p<sub>T</sub> > 30 GeV, \|&eta;\| < 4.7 | `Jet.pTcut`, `Jet.Etacut` |
| jet ID | **Run 2 UL TightLepVeto**, recomputed by `H4LTools::PassJetIDv15` exactly as in CMSSW `PFJetIDSelectionFunctor.h`: `RUN2UL16PUPPI` for 2016, `RUN2ULPUPPI` for 2017/18 | |
| PU jet ID | none (it does not exist for PUPPI) | |
| lepton cleaning | &Delta;R > 0.4 from every **selected** (tight, isolated) lepton | `Jet.LepCleanDR` |
| photon cleaning (&gamma; CR) | &Delta;R > 0.4 from the selected photon | `GammaCR.JetCleanDR` |

The TightLepVeto cuts for \|&eta;\| &le; 2.6 (2017/18), or \|&eta;\| &le; 2.4
(2016), are:

- CEMF < 0.8, CHF > 0, CHM > 0, N<sub>const</sub> > 1
- NEMF < 0.9, MUF < 0.8, NHF < 0.9

The forward regions have their own PUPPI-specific cuts; see the function for the
full table.

The v9 path (AK4 CHS, stored `Jet_jetId` bit 2, PU ID tight below 50 GeV) is still
in the code, behind `nanoAODVersion < 15`.

## b tagging

The tagger is **UParTAK4** (`Jet_btagUParTAK4B`), applied to selected jets with
\|&eta;\| < 2.4. The skim stores counts of jets above each working point
(`HZZ2l2qNu_n{Loose,Medium,Tight}BtagJets`) and does **not** cut on them. The
control plots apply the b-veto with the **Loose** WP. AN-2016/325 used CSVv2
Loose.

The working points come from the BTV `Run2-20XX-UL-NanoAODv15` `btagging.json.gz`:

| year | Loose | Medium | Tight |
| --- | --- | --- | --- |
| 2016preVFP | 0.0387 | 0.1847 | 0.5467 |
| 2016postVFP | 0.0400 | 0.1898 | 0.5538 |
| 2017 | 0.0331 | 0.1776 | 0.5755 |
| 2018 | 0.0308 | 0.1610 | 0.5405 |

## Missing transverse momentum

The analysis uses **PuppiMET** (`PuppiMET_pt`, `PuppiMET_phi`) with the
**MET-&phi; (XY) correction** applied (see
[Corrections](corrections.md#met-phi-correction)). The worker holds only this one
MET. The corrected values are written as `pT_MET` and `phi_MET`. PF `MET_*` is not
used.
