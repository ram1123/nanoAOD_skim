# Event selection

The selection has two layers:

1. **The skim** (`post_proc.py` &rarr; `H4LCppModule` &rarr; `H4LTools`) is loose.
   It defines the regions and stores flags.
2. **The analysis cuts** (MET, p<sub>T</sub><sup>&ell;&ell;</sup> > 55, on-Z, b-veto,
   &Delta;&phi;) are applied downstream, in the control-plot cut files and in
   `higgs_combine/`.

This page describes the skim. The downstream regions are on
[Control plots](control-plots.md#control-regions).

## Order of operations

1. **Preselection**, a ROOT entry list built before any Python runs:
    - `2l2v`: `Sum$(Muon_pt > 20) + Sum$(Electron_pt > 20) >= 2`. The threshold is
      the YAML lepton p<sub>T</sub> cut minus 5 GeV, as a margin for Rochester.
    - `1g`: `Sum$(Photon_pt > 50 && Photon_isScEtaEB) >= 1`.
    - `all`: the OR of the two.

    About 83% of data events never reach Python. The output is identical with
    `--noPreselection`.
2. **Golden JSON** (data).
3. **Trigger:** the event is kept if **any** `TriggerChannels` list fires. Each list
   is also stored as a boolean branch, for example `Triggers_HZZ2l2nu_SR`.
4. **MET filters.**
5. **Object selection** (see [Objects](objects.md)).
6. **Region finding**, below. A &gamma;-CR event is only looked for if no
   dilepton candidate was found.

## Triggers

**Only unprescaled HLT paths** are allowed in the lepton lists. A path missing from
a file's HLT menu, such as the 2016 e&mu; paths in some eras, gets a warning in
`beginFile` and counts as not fired. Check any edit with
`trigger_prescale_check/classify_triggers.py`.

The downstream SR and DY CR require `Triggers_HZZ2l2nu_SR`; the e&mu; CR requires
`Triggers_HZZ2l2nu_EMuCR`.

=== "2016 (pre and post VFP)"

    | list | HLT paths |
    | --- | --- |
    | `Triggers_HZZ2l2nu_SR` | IsoMu24, IsoTkMu24, Mu50, Ele27_WPTight_Gsf, Ele115_CaloIdVT_GsfTrkIdT, Photon175, Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, DoubleEle33_CaloIdL_GsfTrkIdVL_MW |
    | `Triggers_HZZ2l2nu_EMuCR` | Mu8_TrkIsoVVL_Ele23_..._IsoVL and Mu23_TrkIsoVVL_Ele8_..._IsoVL (Run2016B&ndash;G); Mu23_..._Ele12_..._DZ, Mu8_..._Ele23_..._DZ, Mu12_..._Ele23_..._DZ (Run2016F&ndash;H); plus the single-lepton paths of the SR list |

    No MuonEG path is both unprescaled and present over the full year, so the e&mu;
    list ORs era-dependent paths.

=== "2017"

    | list | HLT paths |
    | --- | --- |
    | `Triggers_HZZ2l2nu_SR` | IsoMu27, Mu50, Ele35_WPTight_Gsf, Photon200, Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, DoubleEle33_CaloIdL_MW |
    | `Triggers_HZZ2l2nu_EMuCR` | Mu23_TrkIsoVVL_Ele12_..._IsoVL_DZ, Mu12_TrkIsoVVL_Ele23_..._IsoVL_DZ, Mu8_TrkIsoVVL_Ele23_..._IsoVL_DZ, IsoMu27, Mu50, Ele35_WPTight_Gsf, Photon200 |
    | `Triggers_HZZ2l2nu_GammaCR` | Photon50/75/90/120/165_R9Id90_HE10_IsoM (prescaled), Photon200 |

=== "2018"

    | list | HLT paths |
    | --- | --- |
    | `Triggers_HZZ2l2nu_SR` | IsoMu24, Mu50, Ele32_WPTight_Gsf, Ele115_CaloIdVT_GsfTrkIdT, Photon200, Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, Mu37_TkMu27, Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, DoubleEle25_CaloIdL_MW |
    | `Triggers_HZZ2l2nu_EMuCR` | Mu23_TrkIsoVVL_Ele12_..._IsoVL (+ _DZ), Mu12_TrkIsoVVL_Ele23_..._IsoVL_DZ, Mu8_TrkIsoVVL_Ele23_..._IsoVL_DZ, IsoMu24, Mu50, Ele32_WPTight_Gsf, Ele115_CaloIdVT_GsfTrkIdT, Photon200 |

The YAML files also carry older lists (`Triggers_HZZ4l`, `Triggers_HZZ2l2nu`,
`_SingleLep`, `_DiLep`, `_HighPt`). These still feed the "any list fires" skim
decision, but the analysis does not cut on them.

## ee / &mu;&mu; signal region (`passZZ2l2nuSelection`)

`GetZ1_2l2qOR2l2nu()` and then `ZZSelection_2l2nu()`:

| # | requirement | value (YAML) | cut flow row |
| --- | --- | --- | --- |
| 1 | at least one same-flavour pair of tight leptons with m<sub>&ell;&ell;</sub> in the window. **No charge requirement** (AN &sect;4.3) | 12 < m<sub>&ell;&ell;</sub> < 200 (`HZZ2l2nu.MZLepcut`), `RequireOS: false` | `HZZ2l2qNu_cut2l` |
| 2 | Z&#8321; = the pair with m closest to m<sub>Z</sub> = 91.1876; `HZZ2l2qNu_isELE` from its flavour | | |
| 3 | p<sub>T</sub>(&ell;&#8321;), p<sub>T</sub>(&ell;&#8322;) > 25 | `Leading_Lep_pT`, `SubLeading_Lep_pT` | `cutpTl1l2` |
| 4 | \|&eta;<sub>&ell;</sub>\| < 2.5 | `Lep_eta` | `cutETAl1l2` |
| 5 | \|m<sub>&ell;&ell;</sub> &minus; m<sub>Z</sub>\| < 110 (wide, to keep the &alpha;-method sidebands) | `M_ll_Window` | `cutmZ1Window` |
| 6 | p<sub>T</sub><sup>&ell;&ell;</sup> > 10 | `Pt_ll` | `cutZ1Pt` |
| &ndash; | &Delta;R(&ell;&#8321;, &ell;&#8322;) > 0.3, **flag only** (`HZZ2l2nu_passDeltaRl1l2`) | | info |
| 7 | **extra-lepton veto**: `HZZ2l2nu_nExtraLep == 0` (exactly two tight leptons, no veto lepton) | `VetoLep` | `cutExtraLepVeto` |
| 8 | **soft-muon veto**: `HZZ2l2nu_nExtraSoftMu == 0` | `VetoLep.SoftMu_pT: 3` | `cutSoftMuVeto` |
| 9 | min &Delta;&phi;(jet, MET) &ge; `dPhi_jetMET`. This is **0 in the skim**, so it never rejects. The AN cut of 0.5 is applied downstream | `dPhi_jetMET: 0.0` | `cutdPhiJetMET` |

The opposite-charge result is stored as `HZZ2l2qNu_cutOppositeChargeFlag`.

## e&mu; control region (`passZZ2l2nu_emuCR_Selection`)

`GetZ1_emuCR()` then the same `ZZSelection_2l2nu()` (rows 7&ndash;9). The CR has
its own `HZZemuCR_*` counters.

| requirement | value |
| --- | --- |
| exactly 1 tight electron **and** 1 tight muon | |
| p<sub>T</sub>(&ell;&#8321;), p<sub>T</sub>(&ell;&#8322;) > 25, \|&eta;\| < 2.5 | as for the SR |
| \|m<sub>e&mu;</sub> &minus; m<sub>Z</sub>\| < 160 | `HZZ2l2nu.emuCR.M_ll_Window` |
| p<sub>T</sub><sup>e&mu;</sup> > 25 | `HZZ2l2nu.emuCR.Pt_ll` |
| no charge cut; `HZZ2l2qNu_cutOppositeChargeFlag` stored; `HZZ2l2nu_isEMuCR = 1` | |

## &gamma;+jets control region (`passGammaCR_Selection`) { #gammajets-control-region }

This region is run with `--channels 1g` and follows AN-15-208 &sect;4.2 and
&sect;5.2. It is 2017 only, because only `Input_2017.yml` has a `GammaCR` block;
other years raise `KeyError`.

- Exactly one tight barrel photon (see [Objects](objects.md#photons)),
  and the event has no dilepton candidate.
- The photon **plays the Z**: Z&#8321; = photon four-vector with m = m<sub>Z</sub>.
  `ZZSelection_2l2nu()` then applies the same steps as the dilepton regions. Every
  lepton counts as extra, so any veto lepton rejects the event. The steps are the
  jets (cleaned against the photon), the b-tag counts, the VBF tag, &Delta;&phi; and
  m<sub>T</sub>. `HZZ2l2nu_isGammaCR = 1`.
- **Trigger windows** (`GammaCR.TriggerWindows`): exclusive photon-p<sub>T</sub>
  windows, each tied to one path. `gamma_trigWindowFired` is set when that path
  fired. Data get `gamma_prescaleWeight = FullLumi / EffLumi`, and MC get 1.

| photon p<sub>T</sub> window (GeV) | path | EffLumi (fb<sup>&minus;1</sup>) | weight |
| --- | --- | --- | --- |
| 55 &ndash; 82.5 | HLT_Photon50_R9Id90_HE10_IsoM | 0.30 | 138 |
| 82.5 &ndash; 99 | HLT_Photon75_R9Id90_HE10_IsoM | 1.32 | 31.5 |
| 99 &ndash; 135 | HLT_Photon90_R9Id90_HE10_IsoM | 2.63 | 15.8 |
| 135 &ndash; 181.5 | HLT_Photon120_R9Id90_HE10_IsoM | 7.80 | 5.33 |
| 181.5 &ndash; 230 | HLT_Photon165_R9Id90_HE10_IsoM | 28.34 | 1.47 |
| > 230 | HLT_Photon200 (unprescaled) | 41.54 | 1 |

## Jets, b tags and the VBF tag

These are computed for every region by `CountJetsAndBtags()` and
`ZZSelection_2l2nu()`:

- `HZZ2l2qNu_nJets`: selected jets (p<sub>T</sub> > 30, \|&eta;\| < 4.7, TightLepVeto,
  cleaned).
- `HZZ2l2qNu_n{Loose,Medium,Tight}BtagJets`: UParTAK4, \|&eta;\| < 2.4. No cut.
- `HZZ2l2nu_minDPhi_METAK4`: min \|&Delta;&phi;(jet, MET)\| over the selected jets
  (999 if there are none).
- **VBF tag** `HZZ2l2nu_ifVBF` (AN &sect;4.3, `HZZ2l2nu.VBF`). It requires all of:
    - the **two leading** jets have \|&Delta;&eta;<sub>jj</sub>\| > 4 and
      m<sub>jj</sub> > 500 GeV;
    - the Z (dilepton, e&mu; pair or photon) lies between them in &eta;
      (`ZCentrality: true`);
    - no other selected jet lies between them (`CentralJetVeto: true`).

## Z + MET kinematics

| branch | definition |
| --- | --- |
| `HZZ2l2nu_dPhi_ZMET` | \|&Delta;&phi;(Z&#8321;, MET)\|, stored, no cut in the skim |
| `HZZ2l2nu_ZZpT` | \|**p**<sub>T</sub><sup>&ell;&ell;</sup> + **p**<sub>T</sub><sup>miss</sup>\| |
| `HZZ2l2nu_ZZmT` | AN-2016/325 Eq. 6, with the Z mass on the MET leg (not `TLorentzVector::Mt()`) |

\[
m_T^2 = \left[\sqrt{p_{T,\ell\ell}^2 + m_{\ell\ell}^2} + \sqrt{E_T^{\text{miss}\,2} + m_Z^2}\right]^2 - \left|\vec p_{T,\ell\ell} + \vec E_T^\text{miss}\right|^2
\]

## What the skim does **not** cut on

The following are applied downstream:

- MET > 125 (only an informational `MET > 100` row in the cut flow);
- p<sub>T</sub><sup>&ell;&ell;</sup> > 55;
- \|m<sub>&ell;&ell;</sub> &minus; m<sub>Z</sub>\| < 15;
- the b-veto;
- min &Delta;&phi;(jet, MET) > 0.5 and \|&Delta;&phi;(Z, MET)\| > 0.5;
- &Delta;R(&ell;&#8321;, &ell;&#8322;) > 0.3.
