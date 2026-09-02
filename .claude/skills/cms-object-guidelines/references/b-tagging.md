# b tagging — 2&ell;2&nu; b-Veto &amp; BTV Recommendations

> For the **X/H&rarr;ZZ&rarr;2&ell;2&nu; NanoAOD-tools skim**
> (`post_proc.py` &rarr; `modules/H4LCppModule.py` &rarr; C++ `src/H4LTools.cc`;
> cut values in `config/Input_<year>.yml`).
>
> - Eras: **Run 2 UL 2016 / 2017 / 2018** and some **2022**; NanoAOD **v9** and
>   **v15**. Run-3-2023+ material does not apply.
> - This repo reads **DeepJet** (`Jet_btagDeepFlavB`). L/M/T jet counts are stored
>   as `HZZ2l2qNu_n{Loose,Medium,Tight}BtagJets`; the working points come from
>   `config/Input_<year>.yml` `Jet.deepJet_btag.{Loose,Medium,Tight}` &rarr;
>   `H4LTools::InitializeEvtCut`. `data/btag/*.csv` are legacy DeepCSV inputs.
> - No b-tag scale factor is applied in the nominal skim.
> - Re-check every WP / SF against the current 2&ell;2&nu; note / BTV recommendation.

Responsible POG: **BTV** (B-Tagging &amp; Vertexing POG). Channel strategy:
`references/hzz-2l2nu.md`.

## Role in 2&ell;2&nu; (AN-2016/325 &sect;4.4.5, &sect;5.1)  **[HIG / AN-2016-325]**

b-tagged AK4 jets drive a **b-jet veto** that suppresses the top
(t&#773;t / tW) non-resonant background:

- **signal region**: reject any event containing a b-tagged jet;
- **non-resonant &alpha;-method sidebands**: *require* &ge; 1 b-tagged jet
  (`references/hzz-2l2nu.md` &sect;5, `cms-systematics-statistics`).

Not a signal object; not used for tagging a resonance.

- AN-2016/325 (2016): CSVv2 **loose** `> 0.423`, jets `pT > 30 GeV`,
  `\|&eta;\| < 2.4`. b-veto efficiency uncertainty **2&ndash;4%** on the MC-driven
  processes (signal, WZ, ZZ), from a `b`/`c` mistag-probability downgrade of 2% and
  a light-jet upgrade of 11% (envelope; &sect;7.3, &sect;8.1). The CSVv2 tagger and
  its numeric WP are **superseded** — UL/Run 2+3 use DeepJet.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S0 | **CMS AN-2016/325** &sect;4.3, &sect;4.4.5 (b-veto), &sect;5.1 (&alpha;-method b-tag requirement), &sect;7.3 / &sect;8.1 (b-veto systematic 2&ndash;4%) | 2&ell;2&nu; analysis note &rarr; `references/hzz-2l2nu.md` &sect;2, &sect;8 | 2016/legacy — tagger + WP superseded |
| C1 | This repo: DeepJet WPs | `config/Input_<year>.yml` `Jet.deepJet_btag.{Loose,Medium,Tight}` &rarr; `H4LTools::InitializeEvtCut` (`Input_2018.yml`: 0.0494 / 0.2770 / 0.7264 — 2018 UL DeepJet) | repo |
| C2 | This repo: b-jet counting + veto | `src/H4LTools.cc` (`HZZ2l2qNu_n{Loose,Medium,Tight}BtagJets`; `HZZ2l2nu_cutbtag` in `ZZSelection_2l2nu()`) | repo |
| C3 | Legacy DeepCSV CSVs | `data/btag/*.csv` (`setup.sh` copies them to `$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF/`) | repo |
| P1 | BTV entry point | `btv-wiki.docs.cern.ch/ScaleFactors/` (per-era DeepJet L/M/T WP pages + `btagging.json.gz`) | reference |

Not covered by stored sources &rarr; **Authoritative CMS verification required**:
the DeepJet L/M/T numeric working points for **every era in scope** (only 2018 is
populated in the config), the `btagging.json.gz` version, and — if SFs are added —
the application method (fixed-WP `comb`/`mujets` vs shape/`iterativeFit`), the
per-flavour handling, and the cross-year correlation model.

Classification tags: **[BTV official]**, **[HIG / AN-2016-325]**,
**[Implementation]** (this repo), **[Verify]**, **[Repo divergence]**.

---

## 1. Required context

Exact era; NanoAOD version (v9 vs v15 `Jet_btagDeepFlavB` — the branch is present
in both); which b-tag **working point** the veto uses (loose in the 2016 note);
whether a systematics run needs b-tag SFs.

---

## 2. b-jet selection in this repo (C1, C2)

| Requirement | Value | Note |
|-------------|-------|------|
| AK4 jet pre-selection | see `jets.md` | pT / `\|&eta;\|` / jet ID |
| tagger | **DeepJet** `Jet_btagDeepFlavB` | v9 and v15 |
| working point | `Jet.deepJet_btag.{Loose,Medium,Tight}` from the year's config | `Input_2018.yml` = 0.0494 / 0.2770 / 0.7264 (2018 UL); **other eras: [Verify] / not populated** |
| veto | AN-2016/325: reject the event on a **loose** b-tag | **`HZZ2l2nu_cutbtag` is commented out** in `ZZSelection_2l2nu()` — the veto is not applied — **[Repo divergence]** |

---

## 3. Scale factors  **[Verify]**

**No b-tag SF is applied in the nominal skim.** For a systematics production:

- BTV `btagging.json.gz` (`jsonpog-integration` `POG/BTV/*_UL` for Run 2 v9),
  matching **DeepJet** and the era;
- the SF method (fixed-WP vs shape), per-flavour (`b` / `c` / light) treatment,
  and the systematic decomposition per the BTV recommendation for that tagger;
- the 2016-note-style envelope (b/c mistag &minus;2%, light +11%) only if
  reproducing AN-2016/325 directly — otherwise use the current BTV prescription.

---

## 4. Review checklist

1. Era identified; the config `Jet.deepJet_btag` WPs are populated for that era
   (not just 2018) and match a current btv-wiki page.
2. The b-veto is actually applied — `HZZ2l2nu_cutbtag` gates the event in
   `ZZSelection_2l2nu()` (currently commented out).
3. Veto working point matches the 2&ell;2&nu; note (loose) and the
   &alpha;-method sideband uses the same tagger/WP with `&ge; 1` b-tag.
4. If a systematics run: `btagging.json.gz` matches DeepJet + era; SF method,
   flavour split and correlation model resolved.

---

## 5. Evidence summary

| Item | POG / source | Eras | Established? |
|------|--------------|------|--------------|
| b-veto role, loose WP, 2&ndash;4% systematic | HIG / AN-2016/325 | 2016 | method yes; CSVv2 WP superseded |
| DeepJet as the repo tagger; L/M/T counting | implementation | Run 2 UL / v15 | yes |
| DeepJet L/M/T WPs per era | BTV | all | only 2018 populated — **[Verify]** the rest |
| b-veto applied in `ZZSelection_2l2nu()` | implementation | — | **no** — commented out (**[Repo divergence]**) |
| b-tag SF | BTV | all | **not applied** in the nominal skim |
| SF method + uncertainty decomposition | BTV | all | **Authoritative CMS verification required** |

## Last verified

- AN-2016/325 transcription: skill update.
- Repo cross-check: `config/Input_2018.yml` + `src/H4LTools.cc::ZZSelection_2l2nu()`
  at skill-update time.
- Current BTV recommendation: **not consulted — [Verify]**.
