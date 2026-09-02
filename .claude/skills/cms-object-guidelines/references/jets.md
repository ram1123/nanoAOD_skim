# AK4 Jets — 2&ell;2&nu; Handling &amp; JERC Recommendations

> For the **X/H&rarr;ZZ&rarr;2&ell;2&nu; NanoAOD-tools skim**
> (`post_proc.py` &rarr; `modules/H4LCppModule.py` &rarr; C++ `src/H4LTools.cc`;
> cut values in `config/Input_<year>.yml` `Jet:` &rarr; `H4LTools::InitializeJetcut`).
>
> - Eras: **Run 2 UL 2016 / 2017 / 2018** and some **2022**; NanoAOD **v9** and
>   **v15**. Run-3-2023+ tags/maps do **not** apply.
> - AK4 jets in 2&ell;2&nu; give the **jet-multiplicity category** (=0 / &ge;1 / VBF),
>   the **VBF tag jets** (`m_{jj} > 500`, `\|&Delta;&eta;\| > 4`), and the **b-veto**
>   (`references/b-tagging.md`). `min\|&Delta;&phi;(jet, MET)\| > 0.5` over
>   `pT > 30 GeV` jets is a signal-region cut (`references/hzz-2l2nu.md` &sect;2).
> - v15 has **no** stored `Jet_jetId` / PU-ID; `H4LTools::PassJetIDv15` recomputes
>   jet ID from PF energy fractions (`SetJets` in `modules/H4LCppModule.py`).
> - JEC/JER run **only** under `post_proc.py --WithSyst` via nanoAOD-tools
>   `createJMECorrector`; the nominal skim uses NanoAOD jet pT as-is.
> - AK8 fat jets play no role in 2&ell;2&nu; (the `FatJet` block in
>   `H4LCppModule.analyze()` is commented out on this branch).

Responsible POG: **JERC** (Jet Energy Resolution &amp; Corrections, JME POG).
Channel strategy: `references/hzz-2l2nu.md`.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S1 | Jet Energy Scale Recommendations — CMS JME JERC | `https://cms-jme-jerc.docs.cern.ch/recommendations/jes/` | page snapshot 2026-08-29 |
| S2 | Jet Energy Resolution Recommendations — CMS JME JERC | `https://cms-jme-jerc.docs.cern.ch/recommendations/jer/` | page snapshot 2026-08-29 |
| S3 | Jet Veto Map Recommendations — CMS JME JERC | `https://cms-jme-jerc.docs.cern.ch/recommendations/jet-veto-maps/` | page snapshot 2026-08-29 |
| S0 | **CMS AN-2016/325** &sect;4.3 (jet selection, VBF tag, b-veto, `&Delta;&phi;` cut), &sect;8.1 (JES/JER/unclustered-MET systematics) | `references/hzz-2l2nu.md` &sect;2, &sect;4 | 2016/legacy — **[Verify]** for UL/Run 3 |
| C1 | This repo | `config/Input_<year>.yml` `Jet:` (`pTcut`, `Etacut`, `deepJet_btag`) &rarr; `H4LTools::InitializeJetcut`; `H4LTools::PassJetIDv15`; JEC/JER via nanoAOD-tools `createJMECorrector` under `--WithSyst`; PU-ID SF in `modules/JetSFMaker.py` (`--WithSyst` only) | repo |

Not covered by the stored snapshots &rarr; mark **Authoritative CMS verification
required** if the analysis needs them: full Jet ID working-point fraction cuts, PU
jet ID working points / SFs, the per-source JEC uncertainty split, and the JER
hybrid-method numeric matching thresholds.

Classification tags: **[JERC official]**, **[HIG / AN-2016-325]**,
**[Implementation]** (this repo), **[Verify]**.

---

## 1. Required context

1. **Exact era** — 2016preVFP, 2016postVFP, 2017, 2018, 2022preEE, 2022postEE.
2. **Data or simulation** — L2L3Residual is data-only; JER smearing is MC-only.
3. **NanoAOD campaign** — selects the JEC/JER global tag and jet collection:
   - **v9** (Run 2): AK4 **CHS**; JES `Summer19UL*_V5/V7`. `Jet_jetId` stored.
   - **v15** (Run 2 re-nano): AK4 **PUPPI**; JES `Summer20UL*NanoV15_V1`.
     `Jet_jetId` **not stored** &rarr; recompute (`H4LTools::PassJetIDv15`).
4. **Jet acceptance / ID working point** — analysis choice (`config` `Jet:`),
   not fixed by JERC.

---

## 2. 2&ell;2&nu; jet selection (AN-2016/325 &sect;4.3)  **[HIG / AN-2016-325]**

| Requirement | AN-2016/325 | This repo |
|-------------|-------------|-----------|
| jets counted for the category at | `pT > 30 GeV` | `Jet.pTcut = 30` |
| acceptance | `\|&eta;\| < 4.7` (PF-loose ID) | `Jet.Etacut = 4.7` |
| jet ID | PF-loose (2016 note) | v9: `Jet_jetId`; **v15: recomputed** by `H4LTools::PassJetIDv15` from PF fractions |
| VBF tag | two leading jets `pT > 30`, `\|&Delta;&eta;_{jj}\| > 4`, `m_{jj} > 500 GeV`; + dilepton centrality + central-jet veto | `HZZ2l2nu_VBFIndexJet1/2` picks the highest-`m_{jj}` pair with `\|&Delta;&eta;\| > 4` &amp; `m_{jj} > 500` — **centrality &amp; central-jet veto not implemented** (**[Repo divergence]**, `references/hzz-2l2nu.md` &sect;4) |
| b-veto | reject event with a loose b-tag | `references/b-tagging.md` (currently commented out) |
| `min\|&Delta;&phi;(jet, MET)\| > 0.5` | over jets `pT > 30 GeV` | applied in `ZZSelection_2l2nu()` |

The categorization (=0 / &ge;1 / VBF) and its systematic (Stewart&ndash;Tackmann
jet-bin uncertainty) are in `references/hzz-2l2nu.md` &sect;4 /
`cms-systematics-statistics`.

---

## 3. Correction sequence (JEC/JER — `--WithSyst` only)  **[JERC official]**

```
raw jet pT
  -> JEC:  L1FastJet -> L2Relative -> L3Absolute -> L2L3Residual   (L2L3Residual: DATA only, run-based)
  -> JER smearing                                                   (MC only, after full JEC)
  -> propagate the pT change to Type-I MET   (references/met.md)
```

- JER smearing is applied **after** the JECs; for JER uncertainties apply nominal
  JEC first, then JER up/down (S2). **MC only** — data is never JER-smeared.
- L2L3Residual is **run-based** in the Run 2 (NanoAODv15) and Run 3 JSON files.
- In this repo the whole sequence runs only under `post_proc.py --WithSyst` (via
  nanoAOD-tools `createJMECorrector`); the nominal skim takes NanoAOD jet pT as-is.

### 3.1 JES global tags (S1, in-scope eras)

| Era | v9 (AK4 CHS) | v15 (AK4 PUPPI) |
|-----|--------------|-----------------|
| 2016 preVFP (APV) | `Summer19UL16APV_V7` | `Summer20UL16APVNanoV15_V1` |
| 2016 postVFP | `Summer19UL16_V7` | `Summer20UL16NanoV15_V1` |
| 2017 | `Summer19UL17_V5` | `Summer20UL17NanoV15_V1` |
| 2018 | `Summer19UL18_V5` | `Summer20UL18NanoV15_V1` |
| 2022 preEE / postEE | — | `Summer22_22Sep2023_V*` / `Summer22EE_22Sep2023_V*` — **[Verify]** the exact `_V` against the live JERC page |

### 3.2 JER SF tags (S2, MC)

| Era | JER tag |
|-----|---------|
| 2016 preVFP / postVFP | `Summer20UL16APV_JRV5_MC` / `Summer20UL16_JRV5_MC` |
| 2017 | `Summer19UL17_JRV4_MC` |
| 2018 | `Summer19UL18_JRV3_MC` |
| 2022 preEE / postEE | `Summer22_22Sep2023_JRV2_MC` / `Summer22EE_22Sep2023_JRV2_MC` — **[Verify]** |

Payloads: `correctionlib` JSON from `jsonpog-integration` `POG/JME`. Always
re-check tags against the live `cms-jme-jerc.docs.cern.ch` pages before a
production run.

### 3.3 JES uncertainties  **[Verify]**

Use the JME **regrouped / reduced** source scheme (~11 sources/year: `Absolute`,
`Absolute_<year>`, `BBEC1`, `BBEC1_<year>`, `EC2`, `EC2_<year>`, `HF`, `HF_<year>`,
`RelativeBal`, `RelativeSample_<year>`, `FlavorQCD`) or the full ~25-source list.
Confirm the exact list and the cross-year correlation model against the JME JEC
uncertainty documentation for the campaign. AN-2016/325 (&sect;8.1) recomputes MET,
`M_T`, the jet category and the b-tag for each JES/JER/unclustered-MET &plusmn;1&sigma;.

---

## 4. Jet ID

- Full working-point fraction cuts (loose / tight / tightLepVeto per era):
  **[Verify]** — use the JME Jet ID TWiki for the campaign; never invent the cuts.
- **NanoAODv15**: `Jet_jetId` is not stored; this repo recomputes it in
  `H4LTools::PassJetIDv15(i, isPUPPI)` from the PF energy fractions / multiplicities
  copied by `SetJets` (`chEmEF`, `neEmEF`, `chHEF`, `neHEF`, `muEF`,
  `nConstituents`, `chMultiplicity`, `neMultiplicity`). Verify it against the
  legacy Jet ID TWiki for the era.

---

## 5. Jet Veto Maps (S3)  **[JERC official]**

- Apply the **same map to data and the corresponding MC**.
- **Validation required**: jet (&eta;, &phi;) maps before/after, impact on other
  distributions, documented in the analysis note.
- **Run 2**: veto the **jet** whose axis is in a vetoed region (not the whole
  event). Loose selection: `pT > 15 GeV`; tight jet ID; PU jet ID for CHS jets
  `pT < 50 GeV`; `(chEmEF + neEmEF) < 0.9`; jet not overlapping a PF muon
  (`&Delta;R < 0.2`). Use the strictest data-driven map; UL16 additionally applies
  the `h2hot_mc` map.

| Era | Tag | Map name(s) |
|-----|-----|-------------|
| 2016 (APV &amp; non-APV) | `Summer19UL16_V0` | `h2hot_ul16_plus_hbm2_hbp12_qie11` **and** `h2hot_mc` |
| 2017 | `Summer19UL17_V2` | `h2hot_ul17_plus_hep17_plus_hbpw89` |
| 2018 | `Summer19UL18_V1` | `h2hot_ul18_plus_hem1516_plus_hbp2m1` |

- **2022 (Run 3)**: veto is **mandatory** and **rejects the whole event** if any
  jet passing the minimal selection (`pT > 15 GeV`; `tightLepVeto` ID — recompute
  for v15; `(chEmEF + neEmEF) < 0.9`) is in a vetoed region. Maps
  `Summer22_23Sep2023_V1` (preEE) / `Summer22EE_23Sep2023_V1` (postEE, flagged
  "EE+ leak"). **[Verify]** against the live page.
- This repo does not apply jet veto maps in the nominal skim — **[Verify]** whether
  the 2&ell;2&nu; analysis requires them.

---

## 6. Pileup mitigation / PU jet ID  **[Verify]**

- Run 2 **CHS** (v9): PU jet ID matters for `pT < 50 GeV`. Working points and SFs:
  JME PU jet ID TWiki. This repo applies PU-ID SFs only via `modules/JetSFMaker.py`
  under `--WithSyst`.
- Run 2 / Run 3 **PUPPI** (v15): pileup handled by the PUPPI weights; no separate
  legacy PU jet ID WP.
- The `2.5 < \|&eta;\| < 3.0` data/MC jet-&eta; mismodelling (Run 3) has ad-hoc JME
  mitigations — (a) raise the jet pT threshold there (e.g. `> 50 GeV`); (b) omit
  the stochastic JER component for unmatched jets (needs JME/JERC approval). This
  repo applies **neither**; **[Verify]** whether 2&ell;2&nu; needs one for 2022.

---

## 7. Review checklist

1. Era / data-MC / NanoAOD version identified; jet collection (CHS v9 / PUPPI v15)
   consistent with the JES/JER tag.
2. If `--WithSyst`: JEC level list correct (L2L3Residual data-only, run-based);
   JES tag matches &sect;3.1; JER smearing MC-only after JEC, tag matches &sect;3.2;
   Type-I MET re-propagated (`references/met.md`).
3. v15 jet ID: `H4LTools::PassJetIDv15` checked against the legacy Jet ID TWiki.
4. 2&ell;2&nu; jet cuts: `pT > 30`, VBF `m_{jj} > 500` &amp; `\|&Delta;&eta;\| > 4`
   (+ missing centrality / central-jet veto), `min\|&Delta;&phi;(jet, MET)\| > 0.5`.
5. Jet veto map decision recorded (currently none); if applied, Run 2 = jet veto /
   2022 = event rejection, same map data+MC, before/after validation.
6. JES uncertainty source list + correlation model resolved for a systematics run.

---

## 8. Evidence summary

| Item | POG / source | Eras | Established? |
|------|--------------|------|--------------|
| JEC sequence, run-based L2L3Residual, MC-only hybrid JER | JERC | Run 2 (v9/v15) + 2022 | yes (S1, S2) |
| JES global tags | JERC | 2016&ndash;2018 fixed; 2022 **[Verify]** `_V` | S1 |
| JER SF tags | JERC | 2016&ndash;2018 fixed; 2022 **[Verify]** | S2 |
| v15 jet-ID recompute (`PassJetIDv15`) | implementation | v15 | yes; vs TWiki **[Verify]** |
| Jet veto maps (Run 2 jet veto / 2022 event rejection) | JERC | Run 2 + 2022 | method yes (S3); **not applied** in this skim |
| PU jet ID WPs / SFs | JME | Run 2 CHS | **not in stored sources** — [Verify] |
| JES uncertainty source list | JME | all | **not in stored sources** — [Verify] |
| 2&ell;2&nu; jet selection / VBF tag / `&Delta;&phi;` cut | HIG / AN-2016/325 | 2016 | method yes; numbers **[Verify]** for UL/2022 |

## Last verified

- JERC pages (S1&ndash;S3) snapshotted 2026-08-29 — re-check the live pages before a
  production run.
- AN-2016/325 &sect;4 transcription + repo cross-check (`config/Input_2018.yml`,
  `src/H4LTools.cc`, `modules/H4LCppModule.py`): skill update.
