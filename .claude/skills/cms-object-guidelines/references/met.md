# Missing Transverse Momentum — Stored CMS Recommendations

> This repo is the **H→ZZ→4l / 2l2q / 2l2nu NanoAOD-tools skim**
> (`post_proc.py` → `modules/H4LCppModule.py` → C++ `src/H4LTools.cc`; cut values in
> `config/Input_<year>.yml`).
> - Eras here: **Run 2 UL 2016 / 2017 / 2018** and some **2022**; NanoAOD **v9** and
>   **v15**. Run-3-2022/2023/2024/2025/2026-specific rows below do **not** apply.
> - The CMS-POG recommendation content below is retained as a starting point and
>   **must be re-checked against this analysis's H→ZZ note / HIG group** before it is
>   treated as a requirement.
> - This-repo pointers: **2l2nu is MET-driven** (`HZZ2l2nu_cutMETgT100`, dPhi(jet,MET), `MT_2l2nu`). MET-φ correction **is** applied: `H4LCppModule.analyze()` builds `METPhiCorrector(Campaign.UL_20XX, is_puppi=True)` from nanoAOD-tools and corrects `PuppiMET` → `pT_MET`/`phi_MET` (UL 2016/2017/2018 only). MET filters: `modules/METFilters.py` `passFilters(event, year)`.


Responsible POG: **JME** (JetMET — MET group). Channel strategy: `references/hzz-2l2nu.md`.

MET object used throughout this analysis: **PuppiMET** (`PuppiMET_pt`, `PuppiMET_phi`,
`PuppiMET_sumEt`), with the MET‑φ correction applied in `H4LCppModule.analyze()`.

**The 2l2nu channel is MET‑driven** — MET *is* the Z→νν leg. The signal region is
defined by `MET > 125 GeV` plus `min|Δφ(jet, MET)| > 0.5` and `|Δφ(Z, MET)| > 0.5`,
and the fit observable is the transverse mass `M_T(ℓℓ, MET)` (`references/hzz-2l2nu.md`
§2–§3, AN‑2016/325 §4.3–§4.4, §6.1). The "not MET‑driven / category‑variable only"
description below is **H→μμ (copperhead) legacy** and does not apply to 2l2nu — see §7.
(The 4l channel is genuinely not MET‑driven.)

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S1 | Pointer to the official MET filter list | `docs/Run3_all_basic_Information.md` (MET / Noise Filters) → `https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Run_3_2022_and_2023_data_and_MC` | local review 2026‑08‑31 |
| S2 | **CMS AN‑2016/325** §4.3 (MET flavour + φ correction), §4.4 (Δφ cuts), §6.1 (MET > 125 GeV optimization) | 2l2nu analysis note — transcribed in `references/hzz-2l2nu.md` §2–§3 | 2016/legacy — method only, **[Verify]** numbers |
| C1 | Implementation | `src/copperhead_processor.py` — PuppiMET blocks; `compute_jet_veto_jetfilter` (~L662); PuppiMET–jet horn recipe (~L489) | 2026‑08‑31 |
| C2 | This‑repo 2l2nu MET path | `modules/H4LCppModule.py` `analyze()` `METPhiCorrector(Campaign.UL_20XX, is_puppi=True)` → `pT_MET`/`phi_MET`; `src/H4LTools.cc::ZZSelection_2l2nu()` (`HZZ2l2nu_cutdPhiJetMET`, `HZZ2l2nu_cutMETgT100`, `MT_2l2nu`) | skill‑update |
| — | Type‑I propagation | via the JEC/JER sequence in `jets.md` §2 | 2026‑08‑31 |

Not covered by stored sources → **Authoritative CMS verification required**: the exact
`Flag_*` filter list per era (data vs MC), MET‑φ (XY) correction applicability to
PuppiMET, and any PuppiMET scale/resolution uncertainty prescription.

Classification tags: **[JME official]**, **[Analysis‑specific]**, **[Implementation]**,
**[Verify]**.

---

## 1. Required context

Run 2 vs Run 3; exact era; data vs MC; the MET flavour (**PuppiMET** here, not
`MET`/`CaloMET`); whether the analysis needs MET‑based uncertainties.

---

## 2. Type‑I / propagation

- The AK4 JEC and JER changes applied under `--WithSyst` are propagated to `PuppiMET`
  (`jets.md` §2). **[Implementation]**
- **A MET‑φ (XY) correction IS applied** in this skim for 2l2nu:
  `H4LCppModule.analyze()` builds `METPhiCorrector(Campaign.UL_20XX, is_data=…,
  is_puppi=True)` from nanoAOD‑tools and corrects `PuppiMET.pt/phi` with
  `npv=event.PV_npvs, run=event.run` → `pT_MET` / `phi_MET`. **`Campaign` covers
  UL 2016/2017/2018 only** — any other `year` leaves the corrector undefined and
  the job raises. **[Verify]** whether the current JME recommendation applies a
  φ correction to *PUPPI* MET for the target era (often small/absent for PUPPI —
  confirm, do not assume) and whether the nanoAOD‑tools payload matches it.

---

## 3. MET noise / event filters  **[Verify]**

- The official list to apply is the one linked from S1.
- The exact `Flag_*` branches applied were **not located** in
  `src/copperhead_processor.py` during this review. Confirm which of the recommended
  filters are applied for each era and on data vs MC — typically:
  `goodVertices`, `globalSuperTightHalo2016Filter`,
  `EcalDeadCellTriggerPrimitiveFilter`, `BadPFMuonFilter`, `BadPFMuonDzFilter`,
  `hfNoisyHitsFilter`, `eeBadScFilter` (data only), `ecalBadCalibFilter`.

---

## 4. Analysis‑specific MET treatments (C1)  **[Analysis‑specific / Implementation]**

- **PuppiMET–jet "horn" recipe** (~L489): reject the event if `PuppiMET_pt > 100 GeV`
  **and** there is ≥ 1 AK4 jet with `Δφ(PuppiMET, jet) > 2.9` (plus an EM‑fraction
  condition on the jet). Targets fake MET from the 2.5 < \|η\| < 3.0 jet mismodelling —
  related to the JER known issue in `jets.md` §4.3.
- **Jet‑veto‑map jet filter** (`compute_jet_veto_jetfilter`, ~L662): when the
  jet‑veto‑map *jet* filter removes jets, `PuppiMET_pt` and `PuppiMET_sumEt` are
  recomputed; `PuppiMET_phi` is deliberately **not** reset (avoids a spurious peak at
  φ = 0).

---

## 5. Review checklist

1. MET flavour = PuppiMET everywhere it is used (note the divergence from AN‑2016/325
   PF Type‑I MET — §7).
2. JEC/JER changes propagated to MET (Type‑I) under `--WithSyst`.
3. MET‑φ correction: `METPhiCorrector` `Campaign` matches the era (UL 16/17/18 only)
   and the payload matches the current JME PUPPI‑MET recommendation (§2).
4. Recommended MET filter flags identified for the era and confirmed applied on data
   and MC via `modules/METFilters.py` (§3).
5. **2l2nu (§7)**: `min|Δφ(jet,MET)|>0.5`, `|Δφ(Z,MET)|>0.5`, `MET>125 GeV`
   (rejecting), and the `M_T` definition all present and matching AN‑2016/325.
6. Systematics: lepton‑scale → MET propagation; JES/JER/unclustered‑MET ±1σ
   recompute MET+`M_T`+category+b‑tag; PuppiMET scale/resolution prescription
   resolved if MET‑based uncertainties are needed.

---

## 6. Cross‑check vs this repo's config (as of 2026‑08‑31)

| Observation | Detail |
|-------------|--------|
| MET filters | recommended `Flag_*` list — verify applied per era via `modules/METFilters.py` `passFilters(event, year)` |
| MET‑φ correction | **applied** for 2l2nu in `H4LCppModule.analyze()` (`METPhiCorrector`, UL 2016/17/18 only); a `year` outside that set leaves the corrector undefined and raises |
| MET reset on jet‑veto filter | copperhead‑only; no equivalent in this skim |

---

## 7. 2l2nu MET selection and observable (S2, C2)  **[HIG / AN‑2016‑325]**

For the **2l2nu** channel MET builds the Z→νν leg. Full chain and the divergences
between AN‑2016/325 and the current skim are in `references/hzz-2l2nu.md` §2–§3,
§9; the MET‑specific points:

- **MET flavour**: AN‑2016/325 uses **PF Type‑I MET** (CHS) with the 2016 MET‑φ
  recipe; PUPPI MET was only "stored / explored". **This repo uses PuppiMET +
  `METPhiCorrector`** — an intentional modernization, but it is a divergence from
  the note. **[Verify]** the PuppiMET φ‑correction (`METPhiCorrector` `Campaign`)
  against the current JME PUPPI‑MET recommendation for the target era.
- **Signal‑region cuts** (AN‑2016/325 §4.4, §6.1):
  1. `min|Δφ(jet, MET)| > 0.5` over jets `pT > 30 GeV` — repo `HZZ2l2nu_dPhi_jetMET`
     / `HZZ2l2nu_cutdPhiJetMET`, **applied**.
  2. `|Δφ(Z, MET)| > 0.5` — **not implemented** in `ZZSelection_2l2nu()`.
  3. final `MET > 125 GeV` (optimized, common to all categories) — repo only
     **increments** `HZZ2l2nu_cutMETgT100` at `PuppiMET_pt > 100` and does **not**
     reject the event.
- **Transverse mass** (AN‑2016/325 Eq. 6): the Z→νν leg is given the **PDG Z
  mass**. Repo `MT_2l2nu = (Z1 + Z2_met).Mt()` with `Z2_met` a **massless**
  4‑vector — a different definition. Confirm which is intended before it feeds a
  fit.
- **Systematics** touching MET (AN‑2016/325 §8.1): lepton momentum scale is
  propagated to MET; JES / JER / **unclustered‑MET scale** are varied ±1σ and MET,
  `M_T`, the jet category and the b‑tag are recomputed each time.

---

## 8. Evidence summary

| Item | POG | Eras | Source | Established? |
|------|-----|------|--------|--------------|
| PuppiMET as the MET object (2l2nu) | analysis | UL 16/17/18 | C2 | yes — repo choice; **divergence** from AN‑2016/325 PF Type‑I MET |
| Type‑I propagation of JEC/JER | JME | all | `jets.md` §2 | under `--WithSyst` only |
| MET‑φ (XY) correction on PuppiMET | JME | UL 16/17/18 | C2 | **applied** (`METPhiCorrector`); payload vs current JME rec **[Verify]** |
| 2l2nu SR cuts: `min|Δφ(jet,MET)|>0.5`, `|Δφ(Z,MET)|>0.5`, `MET>125` | HIG | 2016 | S2 | method yes; repo applies only #1 fully — **[Repo divergence]** |
| `M_T(ℓℓ,MET)` = AN‑2016/325 Eq. 6 | HIG | 2016 | S2 | note yes; repo uses a massless‑MET `TLorentzVector::Mt()` — **divergence** |
| MET noise / event filter list | JME | all | S1, `modules/METFilters.py` | verify per era on data/MC |
| PuppiMET scale/resolution uncertainty | JME | all | — | **Authoritative CMS verification required** |

## Last verified

- Local source review: 2026‑08‑31
- AN‑2016/325 (S2) transcription: skill update (this edit)
- Current POG recommendation: pending
