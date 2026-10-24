# Output

Each job writes one hadd-ed skimmed NanoAOD file (`Events`, `Runs`, `LuminosityBlocks`
trees), a `cutFlow` TH1F inside that file, and `cutFlow.json`. Only events that
pass at least one region are kept.

The input branches kept are built at run time from the lists in
`modules/keep_and_drop_list.py`:

- `run`, `event`, `lumi*`, PV, pileup and the physics objects;
- for MC, `LHE*` and `PSWeight*`;
- not `HLT*` or `Flag*`. The trigger decision is stored as the
  `Triggers_HZZ2l2nu_*` booleans.

The `Runs` tree keeps `genEventSumw`, which is needed for normalisation.

## Analysis branches (2&ell;2&nu;)

| group | branches |
| --- | --- |
| region flags | `passZZ2l2nuSelection`, `passZZ2l2nu_emuCR_Selection`, `passGammaCR_Selection`, `HZZ2l2nu_isEMuCR`, `HZZ2l2nu_isGammaCR`, `HZZ2l2qNu_isELE` (Z&#8321; = ee) |
| triggers | one boolean per `TriggerChannels` key, e.g. `Triggers_HZZ2l2nu_SR`, `Triggers_HZZ2l2nu_EMuCR`, `Triggers_HZZ2l2nu_GammaCR` |
| flags (no cut) | `HZZ2l2qNu_cutOppositeChargeFlag`, `HZZ2l2nu_passDeltaRl1l2`, `HZZ2l2nu_ifVBF` |
| leptons | `pTL1/etaL1/phiL1/massL1`, `pTL2/...` (p<sub>T</sub>-ordered, FSR-corrected), `DeltaRl1l2`, `Electron_Fsr_*`, `Muon_Fsr_*` |
| Z&#8321; | `massZ1`, `pTZ1`, `etaZ1`, `phiZ1` |
| MET | `pT_MET`, `phi_MET` (PuppiMET, &phi;-corrected); `pTZ2_met`, `phiZ2_met`, `EneZ2_met` |
| ZZ system | `HZZ2l2nu_ZZmT`, `HZZ2l2nu_ZZpT`, `HZZ2l2nu_dPhi_ZMET` |
| vetoes | `HZZ2l2nu_nExtraLep`, `HZZ2l2nu_nExtraSoftMu`, `HZZ2l2nu_nExtraSoftMuPV` |
| jets | `HZZ2l2qNu_nJets`, `HZZ2l2qNu_n{Loose,Medium,Tight}BtagJets`, `HZZ2l2nu_minDPhi_METAK4`, `pTj1/etaj1/phij1/mj1`, `pTj2/...` |
| VBF | `HZZ2l2nu_VBFIndexJet{1,2}`, `HZZ2l2nu_VBFjet{1,2}_{pT,eta,phi,mass}`, `HZZ2l2nu_VBFdijet_{mass,pT,E}`, `HZZ2l2nu_VBFdEta_jj`, `HZZ2l2nu_VBFdPhi_jj`, `HZZ2l2nu_VBFdR_jj` |
| &gamma; CR | `gamma_pt`, `gamma_eta`, `gamma_phi`, `gamma_trigWindowFired`, `gamma_trigWindowPath`, `gamma_prescaleWeight` |
| weights (MC) | `overallEventWeight`, `overallEventWeight_{puUp,puDown,prefireUp,prefireDown}` (see [Corrections](corrections.md#event-weights)) |

Unfilled values are `-999`, and `minDPhi` is 999 when there are no jets. The 4&ell;
and 2&ell;2q branches (`mass4l`, `D_*`, `massZ2`, `HZZ2l2q_*`, ...) are also
written, but stay at their defaults under `--channels 2l2v`.

## Downstream: `higgs_combine/` { #downstream-higgs_combine }

These are standalone ROOT macros, not built by `scram`:

- `make_shapes.cpp` and `make_shapes_final.cpp` build the shape inputs from the
  skims. Cross sections and categories are hard-coded in the macros.
- `make_datacard.cpp` writes the Combine datacard from those yields.

Run them with ROOT, for example `root -l -b -q higgs_combine/make_shapes.cpp`. The
same normalisation contract applies: weight by `overallEventWeight` and normalise
by &Sigma; `Runs/genEventSumw`. Data from several primary datasets must be
de-duplicated first:

```bash
root -l -b -q 'scripts/analysis/removeDuplicates.C+("inputs.txt","out.root")'
```

Here `inputs.txt` lists one file or pattern per line, in dataset-priority order.
