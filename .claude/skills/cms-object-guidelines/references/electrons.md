# Electrons — 2&ell;2&nu; Selection &amp; EGM Recommendations

> For the **X/H&rarr;ZZ&rarr;2&ell;2&nu; NanoAOD-tools skim**
> (`post_proc.py` &rarr; `modules/H4LCppModule.py` &rarr; C++ `src/H4LTools.cc`;
> cut values in `config/Input_<year>.yml` `Electron:` / `HZZ2l2nu:` &rarr;
> `H4LTools::InitializeElecut` / `InitializeHZZ2l2nuCut`).
>
> - Eras: **Run 2 UL 2016 / 2017 / 2018** and some **2022**; NanoAOD **v9** and
>   **v15**. Run-3-2023+ material does not apply.
> - Electrons build the **Z&rarr;ee leg** — a **primary** selection object.
> - ID flag passed to the worker: v9 `Electron_mvaFall17V2Iso_WP90` &rarr; v15
>   `Electron_mvaIso_WP90` (`modules/H4LCppModule.py` `SetElectrons`); the config
>   `Electron.BDTWP` per-&eta;/per-pT points feed `H4LTools`. FSR recovery:
>   `H4LTools::ElectronFsr`.
> - CMS-POG numbers below are a starting point — re-check against the current
>   2&ell;2&nu; analysis note / HIG group before treating them as a requirement.

Responsible POG: **EGM** (E/Gamma POG). Channel strategy: `references/hzz-2l2nu.md`.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S0 | **CMS AN-2016/325** &sect;4.3 Table 5 (electron ID), &sect;4.3 Eq. 1 (isolation), &sect;4.4 (2&ell;2&nu; pre-selection + 3rd-lepton veto), &sect;7.2 / &sect;8.1 (ID+iso, energy-scale uncertainties) | 2&ell;2&nu; analysis note &rarr; `references/hzz-2l2nu.md` &sect;2 | 2016/legacy — **[Verify]** for UL/Run 3 |
| C1 | Per-year electron cuts | `config/Input_<year>.yml` `Electron:` (`pTcut`, `Etacut`, `Loosedxycut`, `Loosedzcut`, `ttH.WP`, `BDTWP.{Low,Med,High}Eta.{Low,High}PT`, `Isocut`) &rarr; `H4LTools::InitializeElecut` | repo |
| C2 | 2&ell;2&nu; leg cuts | `config/Input_<year>.yml` `HZZ2l2nu:` (`Leading/SubLeading_Lep_pT`, `Lep_eta`, `Pt_ll`, `M_ll_Window`) &rarr; `InitializeHZZ2l2nuCut` | repo |
| C3 | ID flag + collection copy | `modules/H4LCppModule.py` `SetElectrons(...)`; selection logic in `src/H4LTools.cc`; FSR in `H4LTools::ElectronFsr` | repo |
| P1 | EGM entry points | `twiki.cern.ch/twiki/bin/view/CMS/EgammaPOG`; Run 2 UL ID + energy-scale pages; `EgammaIDRecipesRun3` | reference |

Per EGM the Run 3 offline MVA electron ID is `mvaEleID-RunIIIWinter22-iso`
(wp80 / wp90) from CMSSW_126X / NanoV11 — the NanoAOD `Electron_mvaIso_WP90` branch
is its **wp90** point. **[Verify]** for the target era: that Winter22 is still the
current EGM training, the Run 3 barrel&ndash;endcap gap definition, and the
reconstruction + ID + energy scale-and-smearing payload/version.

Not covered by stored sources &rarr; **Authoritative CMS verification required**:
whether HZZ&rarr;2&ell;2&nu; prescribes the plain EGM MVA-iso WP90 or an
analysis-specific BDT working point per era, and the reco/ID/energy-SF prescription.

Classification tags: **[EGM official]**, **[HIG / AN-2016-325]**,
**[Implementation]** (this repo), **[Verify]**, **[Repo divergence]**.

---

## 1. Required context

Era; data vs MC; NanoAOD version (**v9** &rarr; `Electron_mvaFall17V2Iso_WP90`,
`pfRelIso03_all`; **v15** &rarr; `Electron_mvaIso_WP90`, migrated branch order in
`SetElectrons`); intended ID working point; whether electron energy corrections /
SFs are needed (the nominal skim applies **none**).

---

## 2. 2&ell;2&nu; electron selection (AN-2016/325 &sect;4.3&ndash;&sect;4.4)  **[HIG / AN-2016-325]**

Electrons form the Z&rarr;ee candidate. Reference values from AN-2016/325
(2016/legacy — **[Verify]** against the current UL / Run 2+3 note):

| Requirement | AN-2016/325 | This repo |
|-------------|-------------|-----------|
| `pT` | > **25 GeV** | `HZZ2l2nu.Leading/SubLeading_Lep_pT = 25` &check; |
| `\|&eta;\|` | < **2.5**, exclude the gap `1.44 < \|&eta;\| < 1.57` | `Electron.Etacut = 2.5`, `HZZ2l2nu.Lep_eta = 2.5` — gap handling: verify in `src/H4LTools.cc` |
| ID | cut-based **tight** (&sect;4.3 Table 5); &sect;4.4 pre-selection text says "medium" — an internal inconsistency in the note | worker gets the EGM **MVA-iso WP90** flag (`SetElectrons`); config also carries `Electron.BDTWP` (per-&eta;/per-pT HZZ BDT points) + `ttH.WP = 0.8` &rarr; `InitializeElecut`. Confirm which is the applied ID |
| Isolation | tight, **cone R = 0.3** (Eq. 1: `[I_ch + max(I_nh + I_γ − A_eff·&rho;, 0)] / pT`) | `Electron_pfRelIso03_all` (R = 0.3 — matches the note), FSR-subtracted, cut **hard-coded `< 0.15`** in `H4LTools::LeptonSelection()` — the config `Electron.Isocut` (0.15) is passed to `InitializeElecut` but **not used** |
| impact parameter | none layered on the tight WP in the 2016 note | config `Electron.Loosedxycut = 0.045`, `Loosedzcut = 0.2` |
| dilepton | `\|m_{ee} − 91\| < 15 GeV`; leptons **not** charge-required; **reject the event if > 2 lepton candidates** | `HZZ2l2nu.M_ll_Window = 0.0` — window **not applied** — **[Repo divergence]** |
| Z `pT` | `p_T^{ee} > 55 GeV` | `HZZ2l2nu.Pt_ll = 10.0` — **[Repo divergence]** |
| 3rd-lepton veto | loose electron `I_rel < 0.15` & `pT > 10 GeV` (+ loose/soft muon, see `muons.md`) | verify the loose definition in `src/H4LTools.cc` |

AN-2016/325 uncertainties on the ee leg: trigger 2%, ID+iso 2% per 2e event,
electron energy scale 0.6% (barrel) / 1.5% (endcap) propagated to MET
(`references/hzz-2l2nu.md` &sect;8, `cms-systematics-statistics`).

---

## 3. Energy corrections  **[Verify]**

The nominal skim applies **no** electron energy scale or smearing. Since electrons
are a kinematic selection object here, a systematics production needs the **EGM
scale-and-smearing** payload for the era (Run 2 EGM UL; Run 3 `electronSS` JSON),
applied to data and simulation separately, with the selection pT thresholds
re-evaluated on the corrected pT.

---

## 4. Scale factors  **[Verify]**

None applied in the nominal skim. For a systematics production: era-specific EGM
**reconstruction** and **ID(+iso)** scale factors from the EGM `electron.json`
payload, with the working point matching the applied ID
(`mvaIso_WP90` / `mvaFall17V2Iso_WP90` or the HZZ BDT point), plus stat + syst
variations.

---

## 5. Review checklist

1. Era / NanoAOD version identified; correct ID branch
   (`mvaFall17V2Iso_WP90` v9 vs `mvaIso_WP90` v15) and `SetElectrons` argument
   order.
2. 2&ell;2&nu; leg: `pT > 25`, `\|&eta;\| < 2.5` + gap `1.44&ndash;1.57` excluded,
   tight ID, tight iso, `\|m_{ee} − 91\| < 15`, `p_T^{ee} > 55` all applied
   (config currently has `M_ll_Window 0`, `Pt_ll 10` — &sect;2).
3. Which ID is actually enforced — EGM MVA-iso WP90 or the config `BDTWP` points.
4. Dilepton is **not** charge-required; event rejected on a 3rd lepton candidate.
5. If a systematics run: EGM scale-smearing + reco/ID SFs added and validated.

---

## 6. Evidence summary

| Item | POG / source | Eras | Established? |
|------|--------------|------|--------------|
| 2&ell;2&nu; electron selection (pT/&eta;/gap/ID/iso, dilepton window, `p_T^{ee}`, 3rd-lepton veto) | HIG / AN-2016/325 &sect;4 | 2016 | method yes; numbers **[Verify]** for UL/Run 3 |
| ID branch per NanoAOD campaign (`mvaFall17V2Iso_WP90` / `mvaIso_WP90`) | EGM | Run 2 UL / v15 | branch yes; "is this the current EGM/HZZ rec" **[Verify]** |
| Energy scale / smearing | EGM | all | **not applied** in the nominal skim; required for a systematics run |
| Reco + ID scale factors | EGM | all | **not applied** in the nominal skim; required for a systematics run |
| Run 3 MVA ID / gap definition / correctionlib version | EGM | 2022 | **Authoritative CMS verification required** |

## Last verified

- AN-2016/325 &sect;4 / &sect;7 / &sect;8 transcription: skill update.
- Repo cross-check: `config/Input_2018.yml` + `modules/H4LCppModule.py` `SetElectrons`
  at skill-update time.
- Current EGM / HZZ recommendation: **not consulted — [Verify]**.
