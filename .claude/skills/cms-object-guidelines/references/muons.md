# Muons — 2&ell;2&nu; Selection &amp; MUO Recommendations

> For the **X/H&rarr;ZZ&rarr;2&ell;2&nu; NanoAOD-tools skim**
> (`post_proc.py` &rarr; `modules/H4LCppModule.py` &rarr; C++ `src/H4LTools.cc`;
> cut values in `config/Input_<year>.yml` `Muon:` / `HZZ2l2nu:` &rarr;
> `H4LTools::InitializeMucut` / `InitializeHZZ2l2nuCut`).
>
> - Eras: **Run 2 UL 2016 / 2017 / 2018** and some **2022**; NanoAOD **v9** and
>   **v15**. Run-3-2023+ material does not apply.
> - Muons build the **Z&rarr;&mu;&mu; leg** — a **primary** selection object.
> - Momentum correction: nanoAOD-tools `muonScaleRes20XX` (Rochester) &rarr;
>   `Muon_corrected_pt`, added for the year in `post_proc.py`. FSR recovery:
>   `H4LTools::MuonFsr` (from the `FsrPhoton` collection).
> - CMS-POG numbers below are a starting point — re-check against the current
>   2&ell;2&nu; analysis note / HIG group before treating them as a requirement.

Responsible POG: **MUO** (Muon POG). Channel strategy: `references/hzz-2l2nu.md`.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S0 | **CMS AN-2016/325** &sect;4.3 Table 6 (muon "tight arbitration" ID + isolation Eq. 2), &sect;4.4 (2&ell;2&nu; pre-selection, soft/loose 3rd-lepton veto), &sect;7.2 / &sect;8.1 (ID+iso, momentum-scale, soft-veto uncertainties) | 2&ell;2&nu; analysis note &rarr; `references/hzz-2l2nu.md` &sect;2 | 2016/legacy — **[Verify]** for UL/Run 3 |
| C1 | Per-year muon cuts | `config/Input_<year>.yml` `Muon:` (`pTcut`, `Etacut`, `Loose/Tight dxy/dz`, `TightTrackerLayercut`, `TightpTErrorcut`, `HighPtBound`, `Isocut`) &rarr; `H4LTools::InitializeMucut` | repo |
| C2 | 2&ell;2&nu; leg cuts | `config/Input_<year>.yml` `HZZ2l2nu:` (`Leading/SubLeading_Lep_pT`, `Lep_eta`, `Pt_ll`, `M_ll_Window`) &rarr; `InitializeHZZ2l2nuCut` | repo |
| C3 | Trigger lists | `config/Input_<year>.yml` `Triggers_HZZ2l2nu*` (large single-&mu; + double-&mu; + cross OR), evaluated by `modules/Helper.py:PassTrig` | repo |
| C4 | Momentum correction | nanoAOD-tools `muonScaleRes2016/2017/2018` producer &rarr; `Muon_corrected_pt`, in `post_proc.py` | repo |
| C5 | Selection logic | `src/H4LTools.cc` muon selection functions; FSR in `H4LTools::MuonFsr` | repo |
| P1 | MUO POG entry point | `twiki.cern.ch/twiki/bin/view/CMS/MuonPOG#User_Recommendations` (+ the era's UL / Run 3 recommendation page) | reference |

Not covered by stored sources &rarr; **Authoritative CMS verification required**
(consult P1 + the current HZZ note): whether HZZ&rarr;2&ell;2&nu; prescribes **tight**
vs **medium** muon ID for the target era, the numeric ID/iso working point, the
reco/tracking-SF requirement, the Run 3 momentum-calibration prescription, and the
correctionlib payload version.

Classification tags: **[MUO official]**, **[HIG / AN-2016-325]**,
**[Implementation]** (this repo), **[Verify]**, **[Repo divergence]**.

---

## 1. Required context

Era (2016preVFP / 2016postVFP / 2017 / 2018 / 2022); data vs MC; NanoAOD version
(**v9** &rarr; stored `Muon_*` branches; **v15** &rarr; migrated branch names);
intended ID/iso working point; whether muon **or** trigger SFs are needed — the
nominal skim applies **none**.

---

## 2. 2&ell;2&nu; muon selection (AN-2016/325 &sect;4.3&ndash;&sect;4.4)  **[HIG / AN-2016-325]**

Muons form the Z&rarr;&mu;&mu; candidate. Reference values from AN-2016/325
(2016/legacy — **[Verify]** against the current UL / Run 2+3 note):

| Requirement | AN-2016/325 | This repo |
|-------------|-------------|-----------|
| `pT` | > **25 GeV** | `HZZ2l2nu.Leading/SubLeading_Lep_pT = 25` &check; |
| `\|&eta;\|` | < **2.4** | `HZZ2l2nu.Lep_eta = 2.5` (applied to muons too) — **[Repo divergence]** |
| ID | **tight** muon: PF muon, `isGlobal` **or** `isTracker`, tight arbitration, track-fit quality (Table 6) | config `Muon:` `TightTrackerLayercut = 5`, `TightpTErrorcut = 0.3`, `Tight dxy/dz = 0.045 / 0.2`, `HighPtBound = 200` &rarr; `InitializeMucut`; verify the ID logic in `src/H4LTools.cc` matches tight, not medium |
| Isolation | tight PF iso, &Delta;&beta;, **cone R = 0.3** (Eq. 2: `[I_ch + max(I_nh + I_γ − 0.5·I_ch^PU, 0)] / pT`) | `Muon_pfRelIso03_all` (R = 0.3, &Delta;&beta; — matches the note), FSR-subtracted, cut **hard-coded `< 0.2`** in `H4LTools::LeptonSelection()` — the config `Muon.Isocut` (0.2) is passed to `InitializeMucut` but **not used** |
| dilepton | `\|m_{&mu;&mu;} − 91\| < 15 GeV`; leptons **not** required opposite charge; **reject the event if > 2 lepton candidates** | `HZZ2l2nu.M_ll_Window = 0.0` — window **not applied** — **[Repo divergence]** |
| Z `pT` | `p_T^{&mu;&mu;} > 55 GeV` | `HZZ2l2nu.Pt_ll = 10.0` — **[Repo divergence]** |
| 3rd-lepton veto | loose muon `I_rel < 0.2`, **or** soft muon `pT > 3 GeV` (+ loose electron, see `electrons.md`) | verify the loose/soft definitions in `src/H4LTools.cc` |

Event-level: number of good primary vertices > 0; &ge; 1 muon trigger-matched (&sect;3).
Note the 2&ell;2&nu; dilepton is **not** an opposite-charge requirement.

AN-2016/325 uncertainties on the &mu;&mu; leg: trigger 2%, ID+iso 2% per 2&mu;
event, muon momentum scale 1% (propagated to MET), soft-muon-veto efficiency
96&ndash;99% (`references/hzz-2l2nu.md` &sect;8, `cms-systematics-statistics`).

---

## 3. Trigger (C3)

The 2&ell;2&nu; path fires on a **large OR** of single-muon, double-muon and
cross-flavour HLT paths listed under `Triggers_HZZ2l2nu*` in
`config/Input_<year>.yml`, evaluated by `modules/Helper.py:PassTrig`; each trigger
key also becomes a boolean output branch. Single-lepton paths recover double-lepton
inefficiency (AN-2016/325 &sect;4.2).

- **[Verify]** against the current HZZ / MUO recommendation for the era: which
  paths belong in the signal OR, the offline matched-muon pT plateau threshold,
  and whether explicit trigger-object &Delta;R matching is applied (AN-2016/325
  uses the double-muon "reference-trigger" efficiency method, &sect;7.1).

---

## 4. Momentum correction (C4)  **[MUO official]**

- nanoAOD-tools **`muonScaleRes2016 / 2017 / 2018`** producer writes
  `Muon_corrected_pt`; added per year in `post_proc.py` (the "muon scale/resolution
  producer" in the module chain). This is the MUO **Rochester** correction for
  Run 2 UL.
- MC and data receive different operations; **never apply MC smearing to data**.
- **2022 / Run 3**: **[Verify]** which producer / MUO scale-smearing JSON the
  chain should use — not established here.

---

## 5. FSR recovery (C5)  **[Analysis-specific]**

`H4LTools::MuonFsr` selects a recovered photon from the `FsrPhoton` collection
(`config/Input_<year>.yml` `FsrPhoton:` — `pTcut`, `Etacut`, `Isocut`, `dRlcut`,
`dRlOverPtcut`), adds its 4-vector to the muon, and folds its energy into
`pfRelIso04_all` **before** the isolation cut. Applied after the base muon
selection.

---

## 6. Scale factors  **[Verify]**

The **nominal skim applies no muon ID / isolation / trigger scale factor** and no
L1-prefiring weight. If SFs are added for a systematics production:

- MUO `muon_Z.json.gz` (`jsonpog-integration` `POG/MUO/*_UL` for Run 2 v9; CAT
  `metadata/MUO/*` for Run 3), with `NUM_/DEN_` keys matching the **exact**
  selection (tight vs medium ID, PF-iso working point, the HLT path in the OR) for
  the era, plus stat + syst variations;
- reco / tracking SF: **[Verify]** whether the target era's MUO recommendation
  requires a separate one for the ID in use;
- Run 2: L1-prefiring weight (`L1PreFiringWeight.{Nom,Up,Dn}`) — **[Verify]**
  whether it is needed for this final state.

---

## 7. Review checklist

1. Era / data-MC / NanoAOD version identified (v9 vs v15 `Muon_*` branch names).
2. 2&ell;2&nu; leg: `pT > 25`, `\|&eta;\| < 2.4` (**muons**, not 2.5), tight muon ID,
   tight PF iso, `\|m_{&mu;&mu;} − 91\| < 15`, `p_T^{&mu;&mu;} > 55` all applied
   (config currently has `Lep_eta 2.5`, `M_ll_Window 0`, `Pt_ll 10` — &sect;2).
3. Dilepton is **not** charged-required; event rejected on a 3rd lepton candidate.
4. `Muon_corrected_pt` (`muonScaleRes20XX`) used everywhere the muon pT enters;
   MC-only smearing; 2022 producer resolved.
5. FSR photon energy folded into `pfRelIso04_all` before the iso cut.
6. If a systematics run: SF keys match the selection and era; reco-SF /
   L1-prefiring decision recorded.

---

## 8. Evidence summary

| Item | POG / source | Eras | Established? |
|------|--------------|------|--------------|
| 2&ell;2&nu; muon selection (pT/&eta;/ID/iso, dilepton window, `p_T^{&mu;&mu;}`, 3rd-lepton veto) | HIG / AN-2016/325 &sect;4 | 2016 | method yes; numbers **[Verify]** for UL/Run 3 |
| Rochester via `muonScaleRes20XX` &rarr; `Muon_corrected_pt` | MUO | Run 2 UL | yes (implementation); 2022 producer **[Verify]** |
| FSR recovery in `H4LTools::MuonFsr` | analysis | all | yes (implementation) |
| Muon ID/iso/trigger SF | MUO | all | **not applied** in the nominal skim |
| Tight vs medium muon ID recommended for HZZ&rarr;2&ell;2&nu; per era | MUO / HIG | all | **Authoritative CMS verification required** |
| Trigger OR contents + matched-muon plateau per era | MUO / HLT | all | **Authoritative CMS verification required** |
| Run 3 momentum calibration / correctionlib version | MUO | 2022 | **Authoritative CMS verification required** |

## Last verified

- AN-2016/325 &sect;4 / &sect;7 / &sect;8 transcription: skill update.
- Repo cross-check: `config/Input_2018.yml` + `src/H4LTools.cc` at skill-update time.
- Current MUO / HZZ recommendation: **not consulted — [Verify]**.
