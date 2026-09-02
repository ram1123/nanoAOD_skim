# H&rarr;ZZ&rarr;2&ell;2&nu; — Channel Analysis Reference

> ## Scope and provenance — read first
>
> This file documents the **X/H &rarr; ZZ &rarr; 2&ell;2&nu; search strategy** so a review
> can check the skim against it. It is the **channel-level** companion to the
> object files (`met.md`, `electrons.md`, `muons.md`, `b-tagging.md`, `jets.md`).
>
> **Primary source: CMS AN-2016/325** ("Search for a spin-zero high mass resonance
> with the Z(&ell;&ell;)Z(&nu;&nu;) final state in 2016"), 13 TeV, 35.9 fb&#8315;&sup1;,
> CMSSW_8_0_25, MiniAODv2. Transcribed into &sect;2&ndash;&sect;8 below.
>
> - AN-2016/325 is a **legacy 2016 (pre-UltraLegacy) note**. Its *method and
>   selection logic* are the reference; its *numeric working points* are
>   superseded:
>   - b tag: CSVv2 `> 0.423` (2016 loose) &rarr; UL/Run 3 use **DeepJet**
>     (`Jet_btagDeepFlavB`, this repo) with per-era L/M/T WPs (`b-tagging.md`).
>   - electron/muon ID: 2016 cut-based tight / "tight arbitration" &rarr; this repo
>     uses the **HZZ MVA / cut-based-with-IP** selection (`electrons.md`,
>     `muons.md`).
>   - MET: 2016 note uses **PF Type-I MET** with the 2016 MET-&phi; recipe &rarr;
>     this repo uses **PuppiMET** + `METPhiCorrector` (`met.md`).
>   Re-verify every number against the current UL / Run 2+3 2&ell;2&nu; analysis
>   note or the HIG working group before treating it as a requirement.
> - This repo pointers: `H4LTools::GetZ1_2l2qOR2l2nu()` + `ZZSelection_2l2nu()` /
>   `ZZSelection_2l2nu_EMu_CR()` in `src/H4LTools.cc`; cut values in
>   `config/Input_<year>.yml` `HZZ2l2nu:` &rarr; `InitializeHZZ2l2nuCut`; cut-flow
>   bins `dynamicCuts_2l2nu` / `dynamicCuts_2l2nu_emu_CR` in
>   `modules/H4LCppModule.py`.

Responsible: **HIG** (Higgs PAG, HZZ subgroup) for the channel strategy; the object
POGs (EGM / MUO / JME / BTV / LUM) for the pieces — see the object files.

Classification tags: **[HIG / AN-2016-325]** (channel-strategy source),
**[Analysis-specific]**, **[Implementation]** (this repo), **[Verify]** (number
needs a current-era source), **[Repo divergence]** (skim differs from the note).

---

## 1. Required context

Era (2016/2016APV/2017/2018/2022); data or MC; NanoAOD version (v9 / v15);
whether the object under review is on the **signal-region** path or the
**e&mu; control region** path; which **jet-multiplicity category** (=0 / &ge;1 /
VBF). If the review concerns a *number* (WP, threshold, K-factor, uncertainty),
identify whether AN-2016/325 is the only source — if so mark **[Verify]**.

---

## 2. Signal-region event selection (AN-2016/325 &sect;4.3&ndash;&sect;4.4)  **[HIG / AN-2016-325]**

Ordered selection ("pre-selection" in the note):

| # | Requirement | Note value | Repo key | Repo status |
|---|-------------|-----------|----------|-------------|
| 1 | HLT | di-lepton OR single-lepton (single recovers di-lepton inefficiency) | `Triggers_HZZ2l2nu*` in `config/Input_<year>.yml` | present (much larger OR list) |
| 2 | &ge; 2 same-flavour leptons, `pT > 25 GeV`, `\|&eta;\| < 2.5(e) / 2.4(&mu;)`, tight ID, tight iso | &sect;4.3 | `HZZ2l2nu.Leading/SubLeading_Lep_pT = 25`, `Lep_eta` | `Lep_eta = 2.5` for **both** flavours — **[Repo divergence]** (&mu; should be 2.4) |
| 3 | dilepton mass window `\|m_{\ell\ell} - 91\| < 15 GeV` (leptons **not** required opposite charge; reject event if > 2 lepton candidates) | &sect;4.3 | `HZZ2l2nu.M_ll_Window` | config has `M_ll_Window: 0.0` — window **not applied** — **[Repo divergence]** |
| 4 | Z candidate `p_T^{\ell\ell} > 55 GeV` | &sect;4.4.3 | `HZZ2l2nu.Pt_ll` | config has `Pt_ll: 10.0` — **[Repo divergence]** (note wants 55) |
| 5 | 3rd-lepton veto: loose electron `I_rel < 0.15` & `pT > 10`; **or** loose muon `I_rel < 0.2`; **or** soft muon `pT > 3 GeV` | &sect;4.4.4 | `H4LTools` extra-lepton logic | verify the loose/soft definitions match |
| 6 | b-jet veto: **no** b-tagged jet (2016 note: CSVv2 loose `> 0.423`), jets `pT > 30`, `\|&eta;\| < 2.4` | &sect;4.3, &sect;4.4.5 | `HZZ2l2nu_cutbtag`, `deepJet_btag` WPs | **b-veto is commented out** in `ZZSelection_2l2nu()` — **[Repo divergence]** |
| 7 | `min \|&Delta;&phi;(jet, MET)\| > 0.5` over jets `pT > 30 GeV` | &sect;4.4.6 | `HZZ2l2nu_dPhi_jetMET`, `HZZ2l2nu_cutdPhiJetMET` | applied |
| 8 | `\|&Delta;&phi;(Z, MET)\| > 0.5` (removes dilepton recoiling along MET) | &sect;4.4.7 | — | **not applied** — **[Repo divergence]** |
| 9 | final **`MET > 125 GeV`** (optimized, common to all categories, &sect;6.1) | &sect;6 | `HZZ2l2nu_cutMETgT100` | code only **counts** `PuppiMET_pt > 100`; no 125 GeV rejection, event not dropped — **[Repo divergence]** |

The note applies MET filters / sample-cleanup filters upstream (&sect;4.3); this repo
does that in `modules/METFilters.py` `passFilters(event, year)`.

### 2.1 e&mu; control region (non-resonant background)

Same as &sect;2 but the two leptons are **opposite flavour** (e&mu;), still in the Z
mass window; used to predict the non-resonant yield and MT shape in the signal
region (&sect;5.1). Repo: `ZZSelection_2l2nu_EMu_CR()` / `dynamicCuts_2l2nu_emu_CR`.

---

## 3. MET and the transverse mass (AN-2016/325 &sect;4.3)  **[HIG / AN-2016-325]**

- MET flavour in the note: **PF Type-I MET** (CHS), with the 2016 MET-&phi;
  correction applied; PUPPI MET stored but "explored in the near term". Final
  result uses standard PF Type-I MET (&sect;4.4). **This repo uses PuppiMET +
  `METPhiCorrector`** — **[Repo divergence]**; see `met.md`.
- Transverse mass (note Eq. 6), the fit observable:

  ```
  M_T^2 = [ sqrt(pT_ll^2 + m_ll^2) + sqrt(MET^2 + m_Z^2) ]^2
          - [ vec(pT_ll) + vec(MET) ]^2
  ```

  i.e. the Z&rarr;&nu;&nu; leg is assigned the **PDG Z mass**. Repo computes
  `MT_2l2nu = (Z1 + Z2_met).Mt()` with `Z2_met` a **massless** 4-vector
  (`SetPtEtaPhiE(MET,0,MET_phi,MET)`) — this is **not** Eq. 6 — **[Repo
  divergence]**; confirm which MT definition the analysis intends.

---

## 4. Event categorization (AN-2016/325 &sect;4.3)  **[HIG / AN-2016-325]**

Jets counted at `pT > 30 GeV`. Three exclusive categories:

| Category | Definition |
|----------|-----------|
| **VBF** | &ge; 2 jets; two leading jets `pT > 30`, `\|&Delta;&eta;_{jj}\| > 4`, `m_{jj} > 500 GeV`; dilepton centrality between the tag jets; **central-jet veto** (no `pT > 30` jet in `[&eta;_min, &eta;_max]` of the tag jets) |
| **&ge; 1 jet** | fails VBF, `&ge; 1` jet |
| **= 0 jet** | no jet `pT > 30` |

Limits are set per category and combined; ggF and VBF production limits derived
separately (&sect;9.2). Repo: `HZZ2l2nu_ifVBF`, `HZZ2l2nu_VBFIndexJet1/2` select the
highest-`m_{jj}` pair with `|&Delta;&eta;| > 4` & `m_{jj} > 500` — **central-jet veto
and dilepton centrality are not implemented** — **[Repo divergence]**.

---

## 5. Background estimation (AN-2016/325 &sect;5)  **[HIG / AN-2016-325]**

| Background | Method | Key numbers / notes |
|-----------|--------|---------------------|
| **Non-resonant** (t&#773;t, tW, WW, W+jets, Z&rarr;&tau;&tau;, + WZ/WWZ/ZZ&rarr;&tau;&tau; leakage) | different-flavour (&alpha;) method: `N_NR_{ll} = &alpha;_l &middot; N_in_{e&mu;}`, with `&alpha;_l = N_out_{ll} / N_out_{e&mu;}` from `m_{\ell\ell}` sidebands (40&ndash;70, 110&ndash;200 GeV), `MET > 70 GeV`, &ge; 1 b-tag | `&alpha;_ee &asymp; 0.37 &plusmn; 0.01`, `&alpha;_{&mu;&mu;} &asymp; 0.68 &plusmn; 0.01`; MT shape taken from e&mu; **data** (MC for VBF). Cross-checked with the k-method (&sect;5.1.1). |
| **Z+jets** (instrumental MET from hadronic-recoil mismeasurement) | &gamma;+jets data-driven: reweight photon `pT` to the Z `pT` per jet-multiplicity bin; subtract genuine-MET MC; use the reweighted &gamma; sample for the MET and MT templates | genuine-MET fraction &asymp; 25% (40%) after (before) `pT` reweighting at `MET = 125 GeV` (&sect;5.2) |
| **Irreducible** (WZ, ZZ, ZVV) | NLO MC + detector simulation, trusted directly (&sect;5.3) | with the K-factors of &sect;6 |

The `\|&Delta;&phi;(Z, MET)\| < 0.5` region showed instrumental-MET mismodelling
(especially &mu;&mu;) &rarr; the cut in &sect;2 #8 was moved to pre-selection (&sect;5.2).

---

## 6. Theory corrections to MC (AN-2016/325 &sect;2.2, &sect;8.2)  **[HIG / AN-2016-325]**

Applied to the **q&#773;q &rarr; ZZ &rarr; 2&ell;2&nu;** (POWHEG) sample:

- **NLO EWK K-factor**: negative, `~ -4%` at low mass falling to `~ -10%` at
  `m_{ZZ} ~ 600 GeV` and high `pT,Z`; differential in `(s&#770;, t&#770;)` from
  Bierweiler / Gieseke (refs [12,13]); **forced to 1 for `m_{ZZ} < 2 m_Z`**
  (on-shell-ZZ approximation).
- **NNLO QCD K-factor**: taken from the H&rarr;ZZ&rarr;4&ell; 2e2&mu; result
  (Grazzini et al.), `~ 1.1&ndash;1.2` vs `m_{ZZ}`; applicability to 2&ell;2&nu;
  flagged for further study in the note.
- **WZ** also receives an NLO EWK K-factor (virtual + photon-induced, LUXqed PDF).

The EWK-correction **uncertainty** (the recoil-&rho; method, `&rho; &lessgtr; 0.3`)
and the NNLO-QCD-K-factor uncertainty live in
`cms-systematics-statistics/references/systematics.md` §3.5 — not repeated here.

**[Verify]** whether this repo applies any of these K-factors — they are **not**
in `src/H4LTools.cc` or `modules/`. If the downstream `higgs_combine/` macros add
them, record where.

---

> The M_T shape fit, MELA width/interference reweighting, the Eq. 14 μ-model, the
> full systematic-uncertainty list (Combine nuisance names + sizes), and the
> `higgs_combine/` cross-check now live in the **`cms-systematics-statistics`**
> skill. This section is a summary; use that skill for datacard / nuisance / limit
> work.

## 7. Signal model (AN-2016/325 &sect;3, &sect;9.2)  **[HIG / AN-2016-325]**

- ggF + VBF, POWHEG production, JHUGen ZZ&rarr;2&ell;2&nu; decay; masses
  200&ndash;3000 GeV; widths &Gamma; = 5 / 10 / 100 GeV.
- MELA matrix-element **reweighting** models width and the three interference
  terms: H&times;continuum (Int1), h&times;H (Int2), h&times;continuum (Int3).
- Signal-strength / interference model (Eq. 14):

  ```
  &sigma;_gg&rarr;(H&rarr;)VV(&mu;, m_VV) =
      (&mu; - sqrt(&mu;)) &middot; &sigma;_{gg&rarr;H&rarr;VV}
    + sqrt(&mu;)         &middot; &sigma;_{gg&rarr;(H&rarr;)VV}
    + (&mu; - sqrt(&mu;)) &middot; &sigma;_{gg&rarr;VV}
  ```

- Statistical treatment: **shape-based fit on M_T**, per category, CL_s 95% upper
  limits vs mass and width. Note: this is the `higgs_combine/` step, not the skim.

---

## 8. Systematic uncertainties (AN-2016/325 &sect;8, Appendix A)  **[HIG / AN-2016-325]**

Signal + MC-driven backgrounds (WZ, ZZ, ZVV, signal):

| Source | Size (2016 note) |
|--------|------------------|
| Integrated luminosity | **2.6%** (2016) |
| Trigger efficiency | 2% (di-e), 2% (di-&mu;) |
| Lepton ID + isolation | 1% per lepton &rarr; **2% per 2e / 2&mu; event** |
| Lepton momentum scale | &mu; 1%; e 0.6% (barrel) / 1.5% (endcap); propagated to MET |
| JES / JER / unclustered-MET scale | recompute MET, M_T, jet category, b-tag per &plusmn;1&sigma; (Tables 18&ndash;19) |
| b-jet veto efficiency | 2&ndash;4% (b/c downgrade 2%, light upgrade 11% envelope) |
| 3rd-lepton (soft-muon) veto | veto efficiency 96&ndash;99% |
| PDF | NNPDF / PDF4LHC |
| Jet-bin categorization (ggF 0/1/&ge;2) | Stewart&ndash;Tackmann; `&Delta;&sigma;^{0j} = sqrt((&Delta;&sigma;^{&ge;0j})^2 + (&Delta;&sigma;^{&ge;1j})^2)` |
| ZZ / WZ cross section | PDF4LHC + (&mu;_R, &mu;_F) &times; 0.5/2; NLO-EWK uncertainty via &rho; categories (&sect;6) |
| MC statistics | per category (Table 20) |

Data-driven backgrounds:

| Source | Size |
|--------|------|
| Non-resonant (Top/WW) &alpha; method | **15%** total on `&alpha;_e`, `&alpha;_&mu;` (stat + MC-closure bias, &le; 13% at `MET = 125`) |
| Z+jets &gamma;+jets method | 3 components: photon-`pT` reweighting (`< 10%`), genuine-MET subtraction, photon-sample statistics |

---

## 9. Cross-check vs this repo (as of the AN-2016/325 review)

| Item | AN-2016/325 | This repo (`config/`, `src/H4LTools.cc`) |
|------|-------------|------------------------------------------|
| `p_T^{\ell\ell}` cut | > 55 GeV | `Pt_ll: 10.0` — **divergence** |
| `\|m_{\ell\ell} - 91\|` window | < 15 GeV | `M_ll_Window: 0.0` (not applied) — **divergence** |
| Muon `\|&eta;\|` in 2&ell;2&nu; block | < 2.4 | `Lep_eta: 2.5` (both flavours) — **divergence** |
| b-jet veto | reject event with a loose b-tag | commented out in `ZZSelection_2l2nu()` — **divergence** |
| Final MET cut | `> 125 GeV`, drops the event | only counts `PuppiMET_pt > 100`; no rejection — **divergence** |
| `\|&Delta;&phi;(Z, MET)\| > 0.5` | pre-selection cut | not implemented — **divergence** |
| MET flavour | PF Type-I MET (+ 2016 &phi; recipe) | PuppiMET + `METPhiCorrector` — **divergence** (documented, see `met.md`) |
| MT definition | Eq. 6 (Z&rarr;&nu;&nu; leg at `m_Z`) | `(Z1 + massless MET 4-vec).Mt()` — **divergence** |
| b tagger | CSVv2 loose `0.423` | DeepJet `Jet_btagDeepFlavB` L/M/T — expected UL/Run 3 change |
| VBF tag | `m_{jj}>500`, `\|&Delta;&eta;\|>4`, centrality, central-jet veto | `m_{jj}>500` & `\|&Delta;&eta;\|>4` only — centrality / central-jet veto missing |
| K-factors (EWK, NNLO QCD on q&#773;qZZ) | applied | not found in skim — **[Verify]** downstream |

Treat every "divergence" as **analysis-specific inconsistency** to raise with the
analyst, not an "official recommendation violation" — several may be deliberate
choices for the current re-implementation, but none is documented.

---

## 10. Evidence summary

| Item | Source | Established? |
|------|--------|--------------|
| 2&ell;2&nu; selection chain (&sect;2), categorization (&sect;4) | AN-2016/325 &sect;4 | yes — for the 2016 note; **[Verify]** numbers for UL/Run 3 |
| MT observable (Eq. 6) | AN-2016/325 &sect;4.3 | yes — note; repo uses a different form |
| &alpha; / k non-resonant method, `&alpha;_ee`/`&alpha;_{&mu;&mu;}` | AN-2016/325 &sect;5.1 | yes — 2016 values |
| &gamma;+jets Z+jets method | AN-2016/325 &sect;5.2 | yes — method |
| q&#773;qZZ NLO EWK + NNLO QCD K-factors | AN-2016/325 &sect;2.2, refs [12,13,16] | yes — note; application in this repo **[Verify]** |
| Systematics list + 2016 sizes | AN-2016/325 &sect;8 | yes — 2016 sizes; re-derive for the target era |
| b tag WP, lepton ID, MET flavour | AN-2016/325 | **superseded** — use `b-tagging.md` / `electrons.md` / `muons.md` / `met.md` |

## Last verified

- AN-2016/325 transcription: on skill update (this file's creation).
- Current 2&ell;2&nu; UL / Run 2+3 analysis note: **not yet consulted — [Verify]**.
- Repo cross-check (&sect;9): against `config/Input_2018.yml` +
  `src/H4LTools.cc::ZZSelection_2l2nu()` at skill-update time.
