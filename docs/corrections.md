# Corrections and weights

The skim runs the selection **once, at nominal**. Corrections that only reweight
events are applied, or stored with their variations. Corrections that change
which events pass the selection need a re-run of the C++ selection for each
variation, which is not implemented yet.

## Summary

| correction | data | MC | where | status |
| --- | :-: | :-: | --- | --- |
| Golden JSON | &check; | | `PostProcessor(jsonInput=...)`, files in `data/golden_json/` | applied |
| MET noise filters | &check; | &check; | `modules/METFilters.py` | applied |
| Rochester muon scale/resolution | &check; | &check; | `muonScaleRes20XX` &rarr; `Muon_corrected_pt` | nominal applied; up/down produced, not kept |
| MET-&phi; (XY) correction on PuppiMET | &check; | &check; | `METPhiCorrector` in `H4LCppModule` | applied |
| FSR recovery | &check; | &check; | `H4LTools` | applied |
| Jet ID recompute (TightLepVeto) | &check; | &check; | `H4LTools::PassJetIDv15` | applied |
| Pileup reweighting | | &check; | `puWeight_UL20XX` | in `overallEventWeight`, &plusmn;1&sigma; stored |
| L1 prefiring | | &check; | `L1PreFiringWeight_{Nom,Up,Dn}` from NanoAOD | in `overallEventWeight`, &plusmn;1&sigma; stored |
| Generator weight | | &check; | `genWeight` | in `overallEventWeight` |
| &gamma;-CR trigger prescale | &check; | | `gamma_prescaleWeight` | stored; applied in plots |
| LHE scale / PDF / PS weights | | &check; | `LHE*`, `PSWeight*` kept | pass-through, envelope downstream |

## Golden JSON

| year | file |
| --- | --- |
| 2016 (both eras) | `Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt` |
| 2017 | `Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt` |
| 2018 | `Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt` |

## MET noise filters

These follow the [JME Run 2 UL recommendation](https://cms-jme-jmar.docs.cern.ch/recommendations/met/noise_filters/)
and are applied to data **and** MC:

- `Flag_goodVertices`
- `Flag_globalSuperTightHalo2016Filter`
- `Flag_HBHENoiseFilter`, `Flag_HBHENoiseIsoFilter`
- `Flag_EcalDeadCellTriggerPrimitiveFilter`
- `Flag_BadPFMuonFilter`, `Flag_BadPFMuonDzFilter`
- `Flag_hfNoisyHitsFilter` (optional in the recommendation; applied here)
- `Flag_eeBadScFilter`
- `Flag_ecalBadCalibFilter`: **2017 and 2018 only**

`Flag_BadChargedCandidateFilter` is not applied, because JME does not recommend it.

## Muon momentum: Rochester

`muonScaleRes2016pre` / `muonScaleRes2016` / `muonScaleRes2017` / `muonScaleRes2018`
run before the analysis module, chosen by the detected era. Every muon cut, the
Z&#8321; mass and the lepton p<sub>T</sub> use `Muon_corrected_pt`. That is why the
ROOT-level preselection uses (lepton p<sub>T</sub> cut &minus; 5 GeV): the
correction is applied after the preselection.

## MET-&phi; correction { #met-phi-correction }

`METPhiCorrector(campaign, is_data, is_puppi=True)` is built once per job and
called per event as `corrector(PuppiMET.pt, PuppiMET.phi, npv=PV_npvs, run=run)`.
It comes from the `nanoAOD-tools` fork.

| year | campaign |
| --- | --- |
| 2016 data | `UL_2016` (run-dependent) |
| 2016preVFP MC | `UL_2016_APV` |
| 2016postVFP MC | `UL_2016_NONAPV` |
| 2017 | `UL_2017` |
| 2018 | `UL_2018` |

For any other year the uncorrected PuppiMET is used, with a one-time warning. The
corrected (p<sub>T</sub>, &phi;) feed the selection (&Delta;&phi;, m<sub>T</sub>) and
the `pT_MET` / `phi_MET` branches.

## Event weights

MC events carry

\[
\texttt{overallEventWeight} = \texttt{genWeight} \times \texttt{puWeight} \times \texttt{L1PreFiringWeight\_Nom}
\]

plus four variations: `overallEventWeight_puUp/_puDown` (PU &plusmn;1&sigma;) and
`overallEventWeight_prefireUp/_prefireDown`. For data all of them are 1.

- **Pileup:** nanoAOD-tools `puWeight_UL2016/UL2017/UL2018` reweights to the UL
  data profile with a **fixed** MC profile. This producer must run **before**
  `H4LCppModule`. `puAutoWeight_20XX` is pre-UL and depends on how files are split
  across jobs, so it is not used. `puWeight_UL2016` has a single 2016 profile, with
  no preVFP/postVFP split.
- **Prefiring:** this is the combined ECAL (2016&ndash;2017) &times; muon
  (2016&ndash;2018) weight. In 2018 it comes from the muon term, so it is **not**
  &asymp; 1.

!!! warning "Normalisation contract"
    Normalise each MC sample by \(\sigma L / \sum \text{genEventSumw}\), summed over
    the `Runs` tree of **all** skim files. Do not normalise by an event count:
    `genWeight` can be about 100. Do not re-apply `genWeight` or `puWeight` on top
    of `overallEventWeight`.

## &gamma;-CR prescale weight

Every `HLT_Photon*` path below `HLT_Photon200` is prescaled in 2017. Each event
is used only through the path assigned to its photon-p<sub>T</sub> window, and data
are weighted by `gamma_prescaleWeight` = FullLumi / EffLumi (a year-averaged
prescale). See [Event selection &rarr; &gamma; CR](event-selection.md#gammajets-control-region).

## Not applied yet { #not-applied-yet }

| correction | blocker |
| --- | --- |
| Muon / electron reco, ID, iso and trigger SFs | POG payloads not wired in yet. They are weight-only, so they do not need a re-selection |
| UParTAK4 b-tag SF | `btagSFProducer` knows only DeepJet/DeepCSV. Needs the BTV `Run2-20XX-UL-NanoAODv15/btagging.json.gz` |
| Electron energy scale and smearing | needs EGM `electronSS` and a per-variation re-selection |
| JES / JER / unclustered-MET variations | `--WithSyst` adds `createJMECorrector` (AK4PFPuppi for v15). It only **stores** shifted branches, and it shifts PF Type-I MET, not the PuppiMET used here |
| Rochester up/down propagation | needs a per-variation re-selection |
| qqZZ NLO-EWK / NNLO-QCD K-factors | the skim keeps no gen-level m<sub>ZZ</sub>, and `GenVarsProducer` is flagged as not working |
| PU jet ID SF | not applicable on v15: PUPPI jets have no `Jet_puId` |

Luminosity, the &alpha;-method, the instrumental-MET (&gamma;+jets) method and MC
statistics are datacard-level items, handled in `higgs_combine/`.
