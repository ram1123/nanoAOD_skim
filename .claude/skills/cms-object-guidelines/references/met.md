# Missing Transverse Momentum — 2&ell;2&nu; (MET-driven) &amp; JME Recommendations

> For the **X/H&rarr;ZZ&rarr;2&ell;2&nu; NanoAOD-tools skim**
> (`post_proc.py` &rarr; `modules/H4LCppModule.py` &rarr; C++ `src/H4LTools.cc`).
>
> - Eras processed: **Run 2 UL 2016 / 2017 / 2018** and some **2022**; NanoAOD
>   **v9** and **v15**. (The MET-filter code in &sect;3 also carries the full Run 3
>   2022-2026 JME table for completeness.)
> - MET object: **PuppiMET** (`PuppiMET_pt/phi/sumEt`), with the **MET-&phi;
>   correction applied** in `H4LCppModule.analyze()`.
> - **2&ell;2&nu; is MET-driven** — MET *is* the Z&rarr;&nu;&nu; leg. Signal region:
>   `MET > 125 GeV`, `min\|&Delta;&phi;(jet, MET)\| > 0.5`, `\|&Delta;&phi;(Z, MET)\| > 0.5`;
>   fit observable `M_T(&ell;&ell;, MET)` (`references/hzz-2l2nu.md` &sect;2&ndash;&sect;3).
> - MET filters: `modules/METFilters.py` `passFilters(event, year, isMC=…)`
>   (&sect;3, JME table).
> - Re-check numbers against the current 2&ell;2&nu; analysis note / JME group.

Responsible POG: **JME** (JetMET, MET group). Channel strategy:
`references/hzz-2l2nu.md`.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S1 | Official MET (noise / event) filter list — CMS JME | `https://cms-jme-jmar.docs.cern.ch/recommendations/met/noise_filters/` (`#run-2`, `#run-3`) | table transcribed into &sect;3, via user |
| S2 | **CMS AN-2016/325** &sect;4.3 (MET flavour + &phi; correction), &sect;4.4 (`&Delta;&phi;` cuts), &sect;6.1 (`MET > 125 GeV` optimization), &sect;8.1 (JES/JER/unclustered-MET systematics) | 2&ell;2&nu; analysis note &rarr; `references/hzz-2l2nu.md` &sect;2&ndash;&sect;3 | 2016/legacy — method; numbers **[Verify]** |
| C1 | This repo | `modules/H4LCppModule.py` `analyze()` — `METPhiCorrector(Campaign.UL_20XX, is_data=…, is_puppi=True)` &rarr; `pT_MET` / `phi_MET`; `src/H4LTools.cc::ZZSelection_2l2nu()` (`HZZ2l2nu_cutdPhiJetMET`, `HZZ2l2nu_cutMETgT100`, `MT_2l2nu`); `modules/METFilters.py` `passFilters(event, year)` | repo |

Not covered by stored sources &rarr; **Authoritative CMS verification required**:
the MET-&phi; (XY) correction applicability to **PuppiMET** for the target era, the
2025-2026 `ecalBadCalibFilter` handling (JME recommendation under review), and any
PuppiMET scale/resolution uncertainty prescription. The `Flag_*` filter list per
era is now covered by S1 (&sect;3).

Classification tags: **[JME official]**, **[HIG / AN-2016-325]**,
**[Implementation]** (this repo), **[Verify]**, **[Repo divergence]**.

---

## 1. Required context

Exact era; data vs MC; MET flavour (**PuppiMET** here, not `MET` / `CaloMET`);
whether the work needs MET-based systematic variations.

---

## 2. MET-&phi; (XY) correction (C1)

`H4LCppModule.analyze()` builds `METPhiCorrector(Campaign.UL_20XX, is_data=…,
is_puppi=True)` from nanoAOD-tools and corrects `PuppiMET.pt/phi` with
`npv=event.PV_npvs, run=event.run` &rarr; the `pT_MET` / `phi_MET` branches and the
worker's `corr_pt` / `corr_phi`.

- **`Campaign` covers UL 2016 / 2017 / 2018 only** — a `year` outside that set
  leaves the corrector undefined and the job raises.
- **[Verify]** whether the current JME recommendation applies a &phi; correction to
  **PUPPI** MET for the era (often small/absent for PUPPI — confirm, do not
  assume) and whether the nanoAOD-tools payload matches it. In particular there is
  no 2022 `Campaign` — resolve before running 2022.

---

## 3. MET noise / event filters (S1, C1)  **[JME official]**

Applied via `modules/METFilters.py` `passFilters(event, year)` (an AND of all
flags). The **Run 2 UltraLegacy** list, from the JME recommendation table (S1):

| `Flag_*` | 2016 UL | 2017 / 2018 UL | Note |
|----------|:-------:|:--------------:|------|
| `goodVertices` | &check; | &check; | primary vertex |
| `globalSuperTightHalo2016Filter` | &check; | &check; | beam halo |
| `HBHENoiseFilter` | &check; | &check; | |
| `HBHENoiseIsoFilter` | &check; | &check; | |
| `EcalDeadCellTriggerPrimitiveFilter` | &check; | &check; | |
| `BadPFMuonFilter` | &check; | &check; | |
| `BadPFMuonDzFilter` | &check; | &check; | |
| `hfNoisyHitsFilter` | &check; | &check; | listed **Optional** |
| `eeBadScFilter` | &check; | &check; | |
| `ecalBadCalibFilter` | — | &check; | **2017-2018 only** — not in the 2016 UL recommendation |
| `BadChargedCandidateFilter` | &cross; | &cross; | **Not recommended** — do **not** apply |

The **Run 3 (2022-2026)** list (S1 `#run-3`):

| `Flag_*` | 2022 / 2023 | 2024 | 2025 / 2026 | Note |
|----------|:-----------:|:----:|:-----------:|------|
| `goodVertices` | &check; | &check; | &check; | |
| `globalSuperTightHalo2016Filter` | &check; | &check; | &check; | |
| `EcalDeadCellTriggerPrimitiveFilter` | &check; | &check; | &check; | |
| `BadPFMuonFilter` | &check; | &check; | &check; | |
| `BadPFMuonDzFilter` | &check; | &check; | &check; | |
| `hfNoisyHitsFilter` | &check; | &check; | &check; | |
| `eeBadScFilter` | &check; | &check; | &check; | |
| `ecalBadCalibFilter` | *see below* | stored flag | **recommendation under review** — **[Verify]** | |
| `HBHENoiseFilter` / `HBHENoiseIsoFilter` | &cross; | &cross; | &cross; | **not required for Run 3** |

**2022-2023 prompt-reconstruction DATA — ECAL bad crystal** (detector ID
838871812): NanoAOD does not carry the ECAL rec hits to re-run
`Flag_ecalBadCalibFilter`, so **do not use the stored flag** for affected
prompt-reco data. Instead, **for data only**, for **runs 362433&ndash;367144**,
reject the event if `PuppiMET_pt > 100 GeV` **and** &ge; 1 AK4 jet satisfies all
of: `Jet_pt > 50`, `-0.5 < Jet_eta < -0.1`, `-2.1 < Jet_phi < -1.8`,
`(Jet_neEmEF > 0.9 or Jet_chEmEF > 0.9)` — **no** `Jet_jetId` requirement. Not for
simulation; not needed for re-reco datasets; estimated good-data loss `< 0.2%`.
The extra 2022 ECAL crystal noted by JME is handled by the jet veto maps and needs
no filter procedure.

- Applied to **both data and MC** (the S1 tables carry no data-only note, except
  the 2022-2023 prompt-reco special selection which is data-only).
- `modules/METFilters.py` matches these tables as of the skill update:
  - **2016 / 2017 / 2018** — the Run 2 UL set; `ecalBadCalibFilter` for 2017/2018
    only; `BadChargedCandidateFilter` **removed** (was applied — "Not recommended").
  - **2022&ndash;2026** — the Run 3 common list + `Flag_ecalBadCalibFilter` (stored
    flag). The module-level `IS_2022_2023_PROMPT_RECO_DATA` switch (default
    `False`) selects the prompt-reco special selection instead, for data
    (`isMC is False`); the caller passes `isMC=self.isMC`. Set the switch only for
    a 2022/2023 prompt-reco data run.
  - 2025-2026 `ecalBadCalibFilter` handling follows the "stored flag" path —
    **[Verify]** once JME finalizes the recommendation.
- The `Flag_*` branches must exist in the NanoAOD campaign (v9 / v15) — confirm
  with `--DEBUG` on a real input file.

---

## 4. Type-I propagation (`--WithSyst` only)

The AK4 JEC / JER changes applied under `post_proc.py --WithSyst` are propagated to
`PuppiMET` (`jets.md` &sect;3). The nominal skim does not re-derive MET. After any
JES/JER/unclustered-MET variation the MET, `M_T`, jet category and b-tag are
recomputed (AN-2016/325 &sect;8.1 — see `cms-systematics-statistics`).

---

## 5. 2&ell;2&nu; MET selection and observable (S2, C1)  **[HIG / AN-2016-325]**

MET is the Z&rarr;&nu;&nu; leg. Full chain and every divergence between AN-2016/325
and the skim: `references/hzz-2l2nu.md` &sect;2&ndash;&sect;3, &sect;9. MET-specific:

- **MET flavour**: AN-2016/325 uses **PF Type-I MET** (CHS) with the 2016 MET-&phi;
  recipe (PUPPI MET only "stored / explored"). **This repo uses PuppiMET +
  `METPhiCorrector`** — a deliberate modernization, but a divergence from the
  note. **[Verify]** against the current JME PUPPI-MET recommendation.
- **Signal-region cuts** (AN-2016/325 &sect;4.4, &sect;6.1):
  1. `min\|&Delta;&phi;(jet, MET)\| > 0.5` over jets `pT > 30 GeV` — the cut is coded
     (`HZZ2l2nu_cutdPhiJetMET`, `if (minDeltaPhi < HZZ2l2nu_dPhi_jetMET)`) but the
     config threshold `HZZ2l2nu.dPhi_jetMET` is **`0.0`**, so it never fires —
     **[Repo divergence]**.
  2. `\|&Delta;&phi;(Z, MET)\| > 0.5` — **not implemented** in `ZZSelection_2l2nu()`.
  3. final `MET > 125 GeV` (optimized, common to all categories) — the code only
     **increments** `HZZ2l2nu_cutMETgT100` at `PuppiMET_pt > 100` and does **not**
     reject the event.
- **Transverse mass** (AN-2016/325 Eq. 6): the Z&rarr;&nu;&nu; leg is given the PDG
  Z mass. Repo `MT_2l2nu = (Z1 + Z2_met).Mt()` with `Z2_met` a **massless**
  4-vector — a different definition. Confirm which is intended before it feeds a
  fit.

---

## 6. Review checklist

1. MET flavour = PuppiMET everywhere (divergence from AN-2016/325 PF Type-I MET —
   &sect;5).
2. MET-&phi;: `METPhiCorrector` `Campaign` matches the era (UL 16/17/18 only; **no
   2022** — resolve); payload vs the current JME PUPPI-MET recommendation (&sect;2).
3. MET filter `Flag_*` set identified for the era and confirmed applied on data and
   MC via `modules/METFilters.py` (&sect;3).
4. 2&ell;2&nu; SR: `min\|&Delta;&phi;(jet,MET)\| > 0.5`, `\|&Delta;&phi;(Z,MET)\| > 0.5`,
   `MET > 125 GeV` (rejecting), and the `M_T` definition all present and matching
   AN-2016/325 (&sect;5 — currently only #1 fully applied).
5. If `--WithSyst`: JEC/JER changes propagated to MET; unclustered-MET &plusmn;1&sigma;
   recomputes MET + `M_T` + category + b-tag; PuppiMET scale/resolution
   prescription resolved if MET-based uncertainties are needed.

---

## 7. Evidence summary

| Item | POG / source | Eras | Established? |
|------|--------------|------|--------------|
| PuppiMET as the MET object | implementation | UL 16/17/18 | yes — repo choice; **divergence** from AN-2016/325 PF Type-I MET |
| MET-&phi; (XY) correction on PuppiMET | JME | UL 16/17/18 | **applied** (`METPhiCorrector`); no 2022 `Campaign`; payload vs current JME rec **[Verify]** |
| MET noise / event filter list | JME | Run 2 UL + Run 3 | **yes** (S1, &sect;3); `modules/METFilters.py` matches the tables; 2025-2026 `ecalBadCalibFilter` **[Verify]** |
| 2022-2023 prompt-reco ECAL special selection | JME | 2022/2023 prompt data | in `modules/METFilters.py`, gated `IS_2022_2023_PROMPT_RECO_DATA` (default off) |
| Type-I propagation of JEC/JER | JME | all | `--WithSyst` only |
| 2&ell;2&nu; SR cuts (`&Delta;&phi;` &times; 2, `MET > 125`) | HIG / AN-2016/325 | 2016 | method yes; repo: `&Delta;&phi;(jet,MET)` threshold is `0.0`, no `&Delta;&phi;(Z,MET)`, no `MET>125` rejection — **[Repo divergence]** |
| `M_T(&ell;&ell;, MET)` = AN-2016/325 Eq. 6 | HIG / AN-2016/325 | 2016 | note yes; repo uses a massless-MET `TLorentzVector::Mt()` — **divergence** |
| PuppiMET scale/resolution uncertainty | JME | all | **Authoritative CMS verification required** |

## Last verified

- AN-2016/325 (S2) transcription + repo cross-check (`modules/H4LCppModule.py`,
  `src/H4LTools.cc::ZZSelection_2l2nu()`, `modules/METFilters.py`): skill update.
- Current JME MET-filter / PUPPI-MET-&phi; recommendation: **not consulted — [Verify]**.
