# Systematic-uncertainty review — X&rarr;ZZ&rarr;2l2nu skim

_Scope: `HZZ_Analysis_2l2nu_dev` branch, NanoAOD v9 & v15, Run 2 UL
2016preVFP / 2016postVFP / 2017 / 2018. Reference: AN-2016/325 and the
`/cms-object-guidelines` + `/cms-systematics-statistics` skills._
_Last updated 2026-09-02. Reviewed via `/coordinate` (physics-reviewer +
code-reviewer); their corrections are folded in below._

> **2016 is effectively non-functional on this branch.** `config/Input_2016.yml`
> is absent, `post_proc.py` has no NanoAOD-v15 token for 2016, and there is no
> preVFP/postVFP split anywhere (year auto-detect maps both to `year=2016`).
> Everything below that mentions 2016 assumes that gap is closed first.

## 1. Status table

| Systematic | Implemented? | Correct? | Applicable? | Action needed | Relevant file/module |
|---|---|---|---|---|---|
| **Pileup reweight** (`puWeight` &uarr;/&darr;) | **Yes** &mdash; `puWeight_UL20XX` runs before `H4LCppModule`; folded into `overallEventWeight_pu{Up,Down}` | **Yes** (after this pass). Previously used `puAutoWeight_20XX` = pre-UL (ReReco) data profile + per-input-file MC profile (job-split dependent) | Yes (all MC) | Done. Open: no preVFP/postVFP split (needs LUM `puWeights.json.gz`); no 2022 payload | `post_proc.py`, `modules/H4LCppModule.py` |
| **L1 prefiring** (`L1PreFiringWeight_{Nom,Up,Dn}`) | **Yes** &mdash; combined weight folded into `overallEventWeight_prefire{Up,Down}` | Yes as a single combined nuisance | Yes &mdash; ECAL: 2016&ndash;2017; **muon: 2016&ndash;2018** | Optional: carry `L1PreFiringWeight_ECAL_*` and `L1PreFiringWeight_Muon_{Up,Dn,StatUp,StatDn,SystUp,SystDn}` as separate nuisances | `modules/H4LCppModule.py` |
| **Nominal event weight** (`genWeight&middot;puWeight&middot;prefire`) | **Yes** &mdash; `overallEventWeight`; data path = 1.0 | Yes, **if** the downstream contract is honoured (normalise by `genEventSumw`, do not re-apply `genWeight`/`puWeight`) | Yes | Update `higgs_combine/make_shapes*.cpp` to that contract (see &sect;4) | `modules/H4LCppModule.py` |
| **Muon momentum scale/res (Rochester)** | Nominal `Muon_corrected_pt` applied; `&hellip;Up/Down_pt` produced, not kept | Nominal OK; variation propagation absent | Yes | **Deferred** &mdash; per-variation re-selection (m_ll window, p_T^ll, MET, M_T) | `post_proc.py`, `modules/keep_and_drop_list.py` |
| **Electron energy scale/smearing** | No | &mdash; | Yes | **Deferred** &mdash; EGM `electronSS` + re-selection | &mdash; |
| **Muon reco/ID/iso SF** | No | &mdash; | Yes | **Implement now as an event weight** (no re-selection). Blocked only on MUO JSON + verified WP key | `modules/H4LCppModule.py` (new weight factor) |
| **Electron reco/ID/iso SF** | No | &mdash; | Yes | **Implement now as an event weight.** Blocked only on EGM JSON + WP key | same |
| **Trigger SF** (dilepton) | No | &mdash; | Yes | **Implement now as an event weight** | same |
| **JES** (regrouped ~11-source) | Only via `--WithSyst` (stores shifted jet/MET branches; jetType now v9/v15-aware) | No &mdash; shifts PF Type-I MET not PuppiMET; no re-selection | Yes (b-veto, VBF tag, &Delta;&phi;(j,MET), category, M_T) | **Deferred** &mdash; H4LTools per-variation loop + PuppiMET shift | `post_proc.py`, `src/H4LTools.cc` |
| **JER** | Same as JES | No | Yes | **Deferred** &mdash; same | same |
| **Unclustered-MET energy** (`scaleumet`) | No | &mdash; | Yes (MET cut, M_T) | **Deferred** &mdash; PuppiMET `unclustEn` shift + re-selection. NB the skim does not yet reject on `MET>125` (only a counter at `PuppiMET_pt>100`) | &mdash; |
| **MET-&phi; (XY) correction** | Yes (nominal only, UL16/17/18; safe fallback other years) | Correct as implemented | Yes | **[Verify]** whether JME prescribes a &phi; correction for **PUPPI** MET in UL (often negligible) and any PuppiMET scale/res uncertainty | `modules/H4LCppModule.py` |
| **b-tag SF** (DeepJet) | No &mdash; algo string set to `deepjet`, module still disabled | N/A &mdash; b-veto is commented out in `ZZSelection_2l2nu`; also the bundled `data/btag/*.csv` are **DeepCSV** legacy files (no DeepJet payload) | Yes, once the b-veto is enabled | **Implement as an event weight** (no re-selection): enable b-veto (`src/H4LTools.cc`) + supply a DeepJet `btagging.json.gz` | `post_proc.py`, `src/H4LTools.cc` |
| **PU-jet-ID SF** | `JetSFMaker` exists, disabled | N/A | v9/CHS only (v15 PUPPI has no `Jet_puId`) | Enable for v9 only if a PU-ID cut is actually applied. `JetSFMaker.analyze()` reads `jet.puId` unconditionally &rarr; would **raise** on v15, not no-op | `modules/JetSFMaker.py` |
| **Luminosity** | No (correct &mdash; datacard only) | &mdash; | Yes | `lnN` in datacard; Run-2 UL per-year + `lumi_13TeV_*` correlation scheme | `higgs_combine/` |
| **LHE scale (&mu;R/&mu;F, 6-pt)** `LHEScaleWeight` | Pass-through; `Runs` sums preserved | OK as ingredient | Yes | Envelope downstream; drop the two anti-correlated points (0.5,2)/(2,0.5); renormalise with `LHEScaleSumw` | `modules/keep_and_drop_list.py` (`keep LHE*`) |
| **PDF + &alpha;S** `LHEPdfWeight` | Pass-through | OK as ingredient | Yes | Hessian/RMS downstream, renormalise with `LHEPdfSumw`; confirm the branch is present for the POWHEG signal | same |
| **Parton shower** `PSWeight` (ISR/FSR) | Pass-through | OK as ingredient | Yes | ISR/FSR split downstream; check `nPSWeight` (1 vs 4) per sample | same |
| **qqZZ / WZ NLO-EWK K-factor** (+ &rho;-method unc.) | **No** | No | Yes (qq&rarr;ZZ, WZ) | This is a **nominal correction** (shape-changing, &minus;4% to &minus;10% vs m_ZZ), not just a nuisance. Currently **un-applicable**: the skim drops gen ZZ/Z kinematics and `GenVarsProducer` is non-working. Add a gen-level weight in `H4LCppModule`, or keep gen m_ZZ / p_T^Z | `src/H4LTools.cc`, `modules/keep_and_drop_list.py`, `modules/GenVarsProducer.py` |
| **qqZZ / ggZZ NNLO-QCD K-factor** (+ scale unc.) | No | No | Yes | Same as above &mdash; needs gen kinematics in the skim or downstream | same |
| **Non-resonant bkg (different-flavour / &alpha; method)** | No (data-driven) | &mdash; | Yes | Downstream / datacard | `higgs_combine/` |
| **Z+jets instrumental MET (&gamma;+jets method)** | No (data-driven) | &mdash; | Yes | Downstream / datacard | `higgs_combine/` |
| **MC statistical** | No | &mdash; | Yes | Combine `autoMCStats` in the datacard | `higgs_combine/make_datacard.cpp` |

## 2. Files changed

### Pass 1 (initial implementation)

**`modules/H4LCppModule.py`** &mdash; `beginFile()` declares 5 `"F"` branches
(`overallEventWeight`, `overallEventWeight_puUp`, `_puDown`, `_prefireUp`,
`_prefireDown`); `analyze()` computes them MC-only after `keepIt = False`, data
path = 1.0; 5 `fillBranch` calls by the `phi_MET` fill; the stale
`# FIXME: Add weight branch` comment replaced by a NOTE.

**`post_proc.py`** &mdash; PU producer moved before `H4LCppModule`; `--WithSyst`
`NameError` fixed (undefined `muonScaleRes()` / `gammaSF()` removed); b-tag SF
algo `deepcsv`&rarr;`deepjet`; WARNING comment on the JME-corrector limitations.

**`modules/keep_and_drop_list.py`** &mdash; NOTE that `Muon_correctedUp_pt` /
`_correctedDown_pt` are produced but not kept.

### Pass 2 (review corrections)

**`post_proc.py`**
- **PU payload fixed:** `puAutoWeight_{2016,2017,2018}` &rarr;
  `puWeight_UL{2016,2017,2018}`. The old producers targeted the **pre-UltraLegacy
  (ReReco)** data pileup profile and rebuilt the MC profile from each input file
  (`"auto"` mode) &mdash; so the weight depended on how files were split across
  jobs. `puWeight_UL20XX` uses the UL data profile with a fixed MC profile
  (`mcPileupUL20XX.root`). Comment records that no preVFP/postVFP split exists in
  this payload and that 2022 still has none.
- **`--WithSyst` jetType:** `jetType="AK4PFchs"` was hard-coded; now
  `"AK4PFPuppi" if "NanoAODv15" in first_file else "AK4PFchs"` in both the MC and
  data branches (v15 jets are PUPPI).

**`modules/H4LCppModule.py`**
- New helper `_get_weight_branch(event, name, default, warn=)` &mdash;
  nanoAOD-tools raises `RuntimeError("Unknown branch ...")` for an absent branch,
  which `getattr(event, name, default)` does **not** catch. The Pass-1 code would
  therefore have **crashed** (not silently fallen back) on any year with no PU or
  no prefiring producer (e.g. 2022). The helper degrades to `default` with a
  one-time warning instead.
- `analyze()` weight block rewritten to use the helper; warns once if `puWeight`
  is missing, or if `L1PreFiringWeight_Nom` is missing on a Run-2 year.
- `beginFile()` comment expanded into an explicit **downstream normalisation
  contract** (see &sect;4).

**`docs/systematics_review.md`** &mdash; this file; corrected per the two reviews.

Not touched (out of scope / deferred): `src/H4LTools.cc`, `include/H4LTools.h`,
`config/Input_*.yml`.

## 3. How the propagated systematics flow

Per MC event, in `analyze()` (via `_get_weight_branch`, MC only):

```
overallEventWeight             = genWeight * puWeight     * L1PreFiringWeight_Nom
overallEventWeight_puUp        = genWeight * puWeightUp   * L1PreFiringWeight_Nom
overallEventWeight_puDown      = genWeight * puWeightDown * L1PreFiringWeight_Nom
overallEventWeight_prefireUp   = genWeight * puWeight     * L1PreFiringWeight_Up
overallEventWeight_prefireDown = genWeight * puWeight     * L1PreFiringWeight_Dn
```

- **Pileup:** `puWeight_UL20XX` (fixed MC profile &rarr; UL data profile, 69.2 mb
  min-bias, &plusmn;4.6%) writes `puWeight` / `puWeightUp` / `puWeightDown` on
  `event` before `H4LCppModule`. The nominal factor is replaced by its
  &plusmn;1&sigma; value in the two PU-variation branches only.
- **L1 prefiring:** NanoAOD stores the combined ECAL&otimes;muon
  `L1PreFiringWeight_{Nom,Up,Dn}`. The nominal factor is replaced in the two
  prefire-variation branches only. **For 2018 this weight is driven by the muon
  term and is not &asymp; 1** &mdash; `prefireUp/Down` do not collapse to nominal.
- Selection is identical for all five weights &mdash; pure weight systematics, one
  nominal event loop. Downstream fills each histogram once per weight.
- All five branches survive output filtering (`keep *Weight*`).

## 4. Double-counting / duplicated-variation check

- **No double counting inside the skim.** Each variation branch shifts exactly
  one independent ingredient (`puWeight` *or* `L1PreFiringWeight`); `genWeight`
  enters `overallEventWeight` once. PU and prefiring are physically independent.
- **The real risk is at the skim &rarr; `higgs_combine` interface.**
  `make_shapes.cpp` currently weights by a bare `"puWeight"` and normalises by
  `xsec*lumi/(nMC - 2*nNeg)`. Switching to `overallEventWeight` **requires doing
  both at once**:
  1. drop the explicit `puWeight` factor (else PU is squared);
  2. replace the denominator with `genEventSumw` = &Sigma; `genWeight` over all
     generated events, taken from the merged `Runs` tree. For samples whose
     `genWeight` is not &plusmn;1 (POWHEG MiNLO, aMC@NLO), `(nPos - nNeg)` &ne;
     &Sigma; `genWeight` &rarr; normalisation bias.
- `genWeight` / `puWeight` / `L1PreFiringWeight_*` are still written as standalone
  branches (`keep *Weight`, `keep L1PreFiring*`) &mdash; for provenance only.
  `overallEventWeight * genWeight` downstream would double-count. This is stated in
  the `beginFile()` comment.
- `LHEScaleWeight` / `LHEPdfWeight` / `PSWeight` are deliberately **not** folded
  in. When used downstream they must be applied as ratios and renormalised to
  `LHEScaleSumw` / `LHEPdfSumw` (scale/PDF) or `genEventSumw` (PS), else the
  sample cross section shifts.

## 5. Still needs analysis-specific input or external values (not code)

### Implementable now as an event weight (no re-selection) — blocked only on a payload

| Systematic | What is still required |
|---|---|
| Muon reco/ID/iso SF, Electron reco/ID/iso SF | MUO/EGM correctionlib JSON for UL16preVFP / 16postVFP / 17 / 18 + **verified WP keys** matching the applied ID (`Muon` tight, `Electron_mvaIso_WP90`); per-event product of per-lepton SFs, correlated up/down; fold into `overallEventWeight` as a new factor + variation branches. |
| Dilepton trigger SF | Measured efficiency/SF (tag-and-probe, or the AN values &mdash; **unverified for UL**) + up/down; event weight factor. |
| b-tag SF (DeepJet, `effb`) | Re-enable the b-veto in `ZZSelection_2l2nu()` (`src/H4LTools.cc`); supply a DeepJet `btagging.json.gz` (bundled CSVs are DeepCSV); `btagSFProducer` shape/`comb` method with `hf`/`lf`/stat/JES-correlated components as an event weight. |
| qqZZ/WZ NLO-EWK and qqZZ/ggZZ NNLO-QCD K-factors (nominal + unc.) | Gen-level weight vs gen m_ZZ / p_T^Z. **Blocked in the skim first**: stop dropping the gen ZZ kinematics (`keep GenPart*` is commented) or apply the K-factor inside `H4LCppModule`. &rho;-method for the EWK&ndash;QCD correlation. |

### Deferred — needs a per-variation re-run of the C++ selection

| Systematic | What is still required |
|---|---|
| JES (regrouped ~11-source), JER (`resj`) | Re-run `ZZSelection_2l2nu` per shifted jet collection (b-veto, VBF tag, min&#124;&Delta;&phi;(j,MET)&#124;, category, M_T) **and** propagate the shift to **PuppiMET** &mdash; `createJMECorrector` only shifts PF Type-I MET. |
| Unclustered-MET energy (`scaleumet`) | PuppiMET `unclustEn` &plusmn;1&sigma; + recompute the MET cut and M_T. |
| Electron energy scale & smearing | EGM shifts + recompute lepton kinematics, m_ll window, p_T^ll, MET, M_T, category. |
| Muon momentum scale/res (Rochester) | `Muon_correctedUp/Down_pt` must re-drive the selection &mdash; keep the full Muon collection or add an in-worker per-variation loop. |

### Downstream / datacard only

| Systematic | What is still required |
|---|---|
| Luminosity | `lnN` in the datacard &mdash; Run-2 UL per-year + `lumi_13TeV_*` correlation scheme. |
| MC statistical | Combine `autoMCStats`. |
| Non-resonant (different-flavour / &alpha; method) | Transfer-factor stat + closure, on an e&mu;-**data** template. |
| Z+jets instrumental MET (&gamma;+jets) | Normalisation `lnN` + shape + photon-sample stat, on a &gamma;+jets-**data** template. |
| Jet-bin categorisation (Stewart&ndash;Tackmann); interference &mu;-model / K-factor scale envelopes | Built in the fit from the pass-through LHE weights (+ gen kinematics for the K-factor). |
| LHE scale / PDF / &alpha;S / PS envelopes | Ingredients are on the tree; the 6-pt scale envelope, PDF Hessian/RMS and ISR/FSR split are built downstream, renormalised with the `Runs`-tree sums, unphysical scale points removed. |
| 2022 pileup weight | LUM Run-3 `puWeights.json.gz`; `H4LCppModule` currently warns once and leaves `puWeight` out of `overallEventWeight` for 2022. |
| 2016 preVFP/postVFP PU split | LUM `puWeights.json.gz` (`Collisions16_UltraLegacy`, preVFP/postVFP keys) &mdash; the nanoAOD-tools UL2016 ROOT payload has no split. |
