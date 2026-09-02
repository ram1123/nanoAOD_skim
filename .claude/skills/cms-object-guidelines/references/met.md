# Missing Transverse Momentum — 2&ell;2&nu; (MET-driven) &amp; JME Recommendations

> For the **X/H&rarr;ZZ&rarr;2&ell;2&nu; NanoAOD-tools skim**
> (`post_proc.py` &rarr; `modules/H4LCppModule.py` &rarr; C++ `src/H4LTools.cc`).
>
> - Eras: **Run 2 UL 2016 / 2017 / 2018** and some **2022**; NanoAOD **v9** and
>   **v15**. Run-3-2023+ material does not apply.
> - MET object: **PuppiMET** (`PuppiMET_pt/phi/sumEt`), with the **MET-&phi;
>   correction applied** in `H4LCppModule.analyze()`.
> - **2&ell;2&nu; is MET-driven** — MET *is* the Z&rarr;&nu;&nu; leg. Signal region:
>   `MET > 125 GeV`, `min\|&Delta;&phi;(jet, MET)\| > 0.5`, `\|&Delta;&phi;(Z, MET)\| > 0.5`;
>   fit observable `M_T(&ell;&ell;, MET)` (`references/hzz-2l2nu.md` &sect;2&ndash;&sect;3).
> - MET filters: `modules/METFilters.py` `passFilters(event, year)`.
> - Re-check numbers against the current 2&ell;2&nu; analysis note / JME group.

Responsible POG: **JME** (JetMET, MET group). Channel strategy:
`references/hzz-2l2nu.md`.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S1 | Official MET (noise / event) filter list | `twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2` (Run 2 UL + Run 3 sections) | reference |
| S2 | **CMS AN-2016/325** &sect;4.3 (MET flavour + &phi; correction), &sect;4.4 (`&Delta;&phi;` cuts), &sect;6.1 (`MET > 125 GeV` optimization), &sect;8.1 (JES/JER/unclustered-MET systematics) | 2&ell;2&nu; analysis note &rarr; `references/hzz-2l2nu.md` &sect;2&ndash;&sect;3 | 2016/legacy — method; numbers **[Verify]** |
| C1 | This repo | `modules/H4LCppModule.py` `analyze()` — `METPhiCorrector(Campaign.UL_20XX, is_data=…, is_puppi=True)` &rarr; `pT_MET` / `phi_MET`; `src/H4LTools.cc::ZZSelection_2l2nu()` (`HZZ2l2nu_cutdPhiJetMET`, `HZZ2l2nu_cutMETgT100`, `MT_2l2nu`); `modules/METFilters.py` `passFilters(event, year)` | repo |

Not covered by stored sources &rarr; **Authoritative CMS verification required**:
the exact `Flag_*` filter list per era (data vs MC), the MET-&phi; (XY) correction
applicability to **PuppiMET** for the target era, and any PuppiMET
scale/resolution uncertainty prescription.

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

## 3. MET noise / event filters (C1)  **[JME official / Verify]**

Applied via `modules/METFilters.py` `passFilters(event, year)`. Confirm the exact
`Flag_*` branches enabled per era and on data vs MC against S1 — the Run 2 UL /
Run 3 recommended set is typically:
`Flag_goodVertices`, `Flag_globalSuperTightHalo2016Filter`,
`Flag_EcalDeadCellTriggerPrimitiveFilter`, `Flag_BadPFMuonFilter`,
`Flag_BadPFMuonDzFilter`, `Flag_hfNoisyHitsFilter`,
`Flag_eeBadScFilter` (data only), `Flag_ecalBadCalibFilter`.

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
  1. `min\|&Delta;&phi;(jet, MET)\| > 0.5` over jets `pT > 30 GeV` —
     `HZZ2l2nu_dPhi_jetMET` / `HZZ2l2nu_cutdPhiJetMET`, **applied**.
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
| MET noise / event filter list | JME | all | S1 — verify the `Flag_*` set per era, data/MC |
| Type-I propagation of JEC/JER | JME | all | `--WithSyst` only |
| 2&ell;2&nu; SR cuts (`&Delta;&phi;` &times; 2, `MET > 125`) | HIG / AN-2016/325 | 2016 | method yes; repo applies only #1 fully — **[Repo divergence]** |
| `M_T(&ell;&ell;, MET)` = AN-2016/325 Eq. 6 | HIG / AN-2016/325 | 2016 | note yes; repo uses a massless-MET `TLorentzVector::Mt()` — **divergence** |
| PuppiMET scale/resolution uncertainty | JME | all | **Authoritative CMS verification required** |

## Last verified

- AN-2016/325 (S2) transcription + repo cross-check (`modules/H4LCppModule.py`,
  `src/H4LTools.cc::ZZSelection_2l2nu()`, `modules/METFilters.py`): skill update.
- Current JME MET-filter / PUPPI-MET-&phi; recommendation: **not consulted — [Verify]**.
