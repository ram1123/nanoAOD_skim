# Electrons — Stored CMS Recommendations

> This repo is the **H→ZZ→4l / 2l2q / 2l2nu NanoAOD-tools skim**
> (`post_proc.py` → `modules/H4LCppModule.py` → C++ `src/H4LTools.cc`; cut values in
> `config/Input_<year>.yml`). 
> - Eras here: **Run 2 UL 2016 / 2017 / 2018**; NanoAOD **v9** and
>   **v15**. Run-3-2022/2023/2024/2025/2026-specific rows below do **not** apply.
> - The CMS-POG recommendation content below is retained as a starting point and
>   **must be re-checked against this analysis's H→ZZ note / HIG group** before it is
>   treated as a requirement.
> - This-repo pointers: cuts in `config/Input_<year>.yml` `Electron:` (incl. `BDTWP` MVA points) → `H4LTools::InitializeElecut`; ID branch is v9 `Electron_mvaFall17V2Iso_WP90` / v15 `Electron_mvaIso_WP90` (`modules/H4LCppModule.py` `SetElectrons`); FSR recovery in `H4LTools::ElectronFsr`. **Primary** selection object (4e / 2e2μ, and Z→ee in 2l2q/2l2nu) — NOT a veto object here.


Responsible POG: **EGM** (E/Gamma POG).

Role in this analysis (**H→ZZ**): electrons are a **primary selection object** — they
build the Z→ee leg in the 4e / 2e2µ (4l), 2l2q and 2l2nu channels. They are **not** a
veto object here. The veto-definition text below is from H→µµ and does not apply;
use `config/Input_<year>.yml` `Electron:` + `H4LTools::InitializeElecut` /
`goodLooseElectrons2012` / `passTight_BDT_Id` as the actual selection.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S1 | This analysis's electron‑veto definition | `docs/Official_recommendation.md` (Electron Selection table) | local review 2026‑08‑31 |
| S0 | **CMS AN‑2016/325** §4.3 Table 5 (electron ID), §4.3 Eq. 1 (isolation), §4.4 (2l2nu pre‑selection, 3rd‑lepton veto) | 2l2nu analysis note → `references/hzz-2l2nu.md` §2, §2b below | 2016/legacy — **[Verify]** for UL/Run 3 |
| S2 | EGM Run 3 offline ID names + WPs (`EgammaIDRecipesRun3`) | CMS PdmV `PdmVRun3Analysis`, "Notes from POGs / From E/Gamma" — see `lumi.md` §9 | via user, 2026‑09‑01 |
| C1 | Per‑year working‑point keys | `configs/parameters/electron.yaml` | 2026‑08‑31 |
| C2 | Implementation | `src/copperhead_processor.py` (electron‑veto block) | 2026‑08‑31 |

Per S2, the Run 3 offline MVA electron ID is `mvaEleID-RunIIIWinter22-iso` (wp80,
wp90), available from CMSSW_126X / NanoV11 — the NanoAOD `Electron_mvaIso_WP90` branch
is its **wp90** point. Still **[Verify]** for the target era: that Winter22 is the
current EGM training (not superseded for 2023/2024/2025), the Run 3 gap definition, the
correctionlib payload/version, and reconstruction + ID scale factors (needed only if
electrons become a selection object).

Classification tags: **[EGM official]**, **[Analysis‑specific]**, **[Implementation]**,
**[Verify]**.

---

## 1. Required context

Run 2 vs Run 3; exact era; NanoAOD campaign (v9 → `mvaFall17V2Iso_WP90`; v12/v15 →
`mvaIso_WP90`); whether electrons stay a veto object or are promoted to a selection
object (changes everything below).

---

## 2. Electron pre‑selection (veto definition)  **[Analysis‑specific]**

An event is **rejected** if any electron passes all of:

| Requirement | Value | Classification |
|-------------|-------|----------------|
| pT | > **20 GeV** | analysis‑specific (`electron_pt_cut`) |
| \|η\| | < **2.5** | detector acceptance (`electron_eta_cut`) |
| barrel–endcap gap | exclude **1.44 < \|η\| < 1.57** | standard EGM exclusion (S1) |
| ID | MVA ID **with isolation**, **WP90** | EGM POG working point |

Working‑point key by NanoAOD campaign (C1):

| NanoAOD | Key | EGM ID (S2) |
|---------|-----|-------------|
| v9 (Run 2 UL) | `mvaFall17V2Iso_WP90` | Fall17V2 MVA, iso, 90 % WP |
| v12 / v15 | `mvaIso_WP90` | `mvaEleID-RunIIIWinter22-iso` **wp90** |

No IP or extra isolation cut is layered on top — the MVA‑with‑iso WP is the whole
definition.

---

## 2b. 2l2nu / 2l2q — electrons as the Z→ee leg (AN‑2016/325 §4.3–§4.4)  **[HIG / AN‑2016‑325]**

For the **2l2nu** channel (the current focus) electrons build the Z→ee candidate;
they are a **primary selection object**, not a veto. Reference numbers from
AN‑2016/325 (2016/legacy — **[Verify]** vs the current UL/Run 2+3 note; see
`references/hzz-2l2nu.md` §2):

| Requirement | AN‑2016/325 value | This repo |
|-------------|-------------------|-----------|
| `pT` | > **25 GeV** | `HZZ2l2nu.Leading/SubLeading_Lep_pT = 25` ✓ |
| `\|η\|` | < **2.5**, exclude gap `1.44 < \|η\| < 1.57` | `HZZ2l2nu.Lep_eta = 2.5` ✓ (gap handling: verify) |
| ID | cut‑based **tight** (§4.3 Table 5); §4.4 text says "medium" — note internal inconsistency | repo uses the **HZZ MVA WP** (`BDTWP` in `config`) / v9 `mvaFall17V2Iso_WP90` vs v15 `mvaIso_WP90` — an intentional modernization |
| Isolation | tight, `I_rel` (Eq. 1: `[I_ch + max(I_nh + I_γ − A_eff·ρ, 0)]/pT`) | `Electron.Isocut = 0.15` |
| dilepton | `\|m_{ee} − 91\| < 15 GeV`, **not** charge‑required, reject if > 2 e candidates | config `M_ll_Window: 0.0` — **not applied** — **[Repo divergence]** |
| Z `pT` | `p_T^{ee} > 55 GeV` | config `Pt_ll: 10.0` — **[Repo divergence]** |
| 3rd‑lepton veto | loose electron `I_rel < 0.15` & `pT > 10 GeV` | verify the loose definition |

Because electrons are now a selection object, **energy scale + smearing and
reco/ID(+iso) scale factors are required** (§3, §4) — not optional as the
veto‑only text implies.

---

## 3. Energy corrections  **[Verify]**

Electron energy scale and smearing (EGM) are **not applied** — electrons only enter a
pass/fail veto. If electrons are ever promoted to a kinematic object:

- apply the EGM scale‑and‑smearing payload for the era (Run 2 EGM UL; Run 3 `electronSS`
  JSON) to data and simulation separately;
- re‑evaluate the veto pT threshold against the corrected pT.

---

## 4. Scale factors  **[Verify]**

None applied (veto‑only). If promoted to a selection object, era‑specific EGM
reconstruction and ID(+iso) scale factors from the EGM `electron.json` payload are
required, with the WP matching `mvaIso_WP90` / `mvaFall17V2Iso_WP90`, plus stat + syst
variations.

---

## 5. Review checklist

1. Era / NanoAOD version identified; correct WP key (`mvaFall17V2Iso_WP90` vs
   `mvaIso_WP90`).
2. pT (20), \|η\| (2.5), gap exclusion (1.44–1.57) applied per C1/C2.
3. Electrons used only as a veto — no energy correction / SF expected; if that changed,
   §3 and §4 now apply.
4. If promoted: EGM scale‑smearing + reco/ID SFs added and validated.

---

## 6. Evidence summary

| Item | POG | Eras | Source | Established? |
|------|-----|------|--------|--------------|
| Veto pre‑selection (pT/η/gap/WP90) | analysis | all | S1, C1 | yes — analysis choice |
| WP key per NanoAOD campaign | EGM | all | C1 | key yes; "is this the current EGM rec" **[Verify]** |
| Energy scale/smearing | EGM | all | — | not applied (veto‑only); **required if promoted** |
| Reco + ID scale factors | EGM | all | — | not applied (veto‑only); **required if promoted** |
| Run 3 MVA ID recommendation / gap definition / correctionlib version | EGM | Run 3 | — | **Authoritative CMS verification required** |

## Last verified

- Local source review: 2026‑08‑31
- Current POG recommendation: pending
