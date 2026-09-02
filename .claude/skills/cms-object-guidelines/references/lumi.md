# Luminosity — 2&ell;2&nu; &amp; LUM Recommendations

> For the **X/H&rarr;ZZ&rarr;2&ell;2&nu; NanoAOD-tools skim**
> (`post_proc.py` &rarr; `modules/H4LCppModule.py` &rarr; C++ `src/H4LTools.cc`).
>
> - Eras: **Run 2 UL 2016 / 2017 / 2018** (13 TeV) and some **2022** (13.6 TeV);
>   NanoAOD **v9** and **v15**. Run-3-2023+ material does not apply.
> - Golden-JSON lumimask applied to **data** by `PostProcessor(jsonInput=...)` —
>   the file is chosen per year in `post_proc.py` from `data/golden_json/`.
> - Per-year integrated luminosity: `config/Input_<year>.yml` `lumi:` (scales MC).
> - `scripts/analysis/dump_lumi.py` + `LumiDumper.py` (a `Module`) collect the
>   processed run &rarr; lumi-section sets and write JSON.
> - Use the **Run 2 13 TeV UltraLegacy** LUM recommendation; re-check any 2022
>   value against the current LUM Run 3 page.

Responsible POG: **LUM** (integrated-luminosity values, their uncertainty, and the
cross-year / datacard-nuisance scheme); **DQM-DC** for the certified golden JSONs.
Channel strategy: `references/hzz-2l2nu.md`.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| S0 | **CMS AN-2016/325** &sect;8.1 — 2016 luminosity uncertainty **2.6%** | 2&ell;2&nu; analysis note &rarr; `references/hzz-2l2nu.md` &sect;8 | 2016/legacy |
| P1 | LUM Run 2 recommendation | `twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2` (UltraLegacy values, per-year uncertainties, `lumi_13TeV_*` correlation scheme) | reference |
| P2 | brilcalc setup | `twiki.cern.ch/twiki/bin/view/CMS/BrilcalcQuickStart` | reference |
| P3 | Certified golden JSONs | `https://cms-service-dqmdc.web.cern.ch/CAF/certification/` | reference |
| C1 | This repo | `post_proc.py` (`jsonInput` = a `data/golden_json/` file per year); `config/Input_<year>.yml` `lumi:`; `scripts/analysis/{dump_lumi.py,LumiDumper.py}` | repo |

Not covered by stored sources &rarr; **Authoritative CMS verification required**:
the exact per-year UL uncertainty and its `lumi_13TeV_*` correlation model (use
P1), and the 2022 luminosity value + uncertainty + golden-JSON file (use the LUM
Run 3 page).

Classification tags: **[LUM official]**, **[HIG / AN-2016-325]**,
**[Implementation]** (this repo), **[Verify]**.

---

## 1. Required context

Exact era; **data vs MC** — the golden JSON is applied to data, the
integrated-luminosity value scales MC; whether the result combines years (&rarr; the
LUM Run 2 correlation scheme, not a flat number).

---

## 2. Golden JSON / lumimask (C1)

Applied to **data only**, as `PostProcessor(..., jsonInput=<file>)`. Standard UL
certifications:

| Era | Golden JSON |
|-----|-------------|
| 2016 (preVFP &amp; postVFP) | `Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt` |
| 2017 | `Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt` |
| 2018 | `Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt` |
| 2022 | `Cert_Collisions2022_355100_362760_Golden.json` — **[Verify]** the file the analysis should use |

Confirm the file `post_proc.py` picks per year from `data/golden_json/` matches
these. Use the **Golden** JSON (not a muon-only certification).

---

## 3. Integrated luminosity for MC normalisation (C1)

The per-year value in `config/Input_<year>.yml` `lumi:` scales the simulated
samples. UL 13 TeV values (fb&#8315;&sup1;, consistent with the LUM UltraLegacy
recommendation):

| Era | fb&#8315;&sup1; |
|-----|-----|
| 2016 preVFP | 19.5 |
| 2016 postVFP | 16.81 |
| 2016 total | 36.31 |
| 2017 | 41.48 |
| 2018 | 59.83 |
| Run 2 total | &asymp; 138 |
| 2022 | **[Verify]** (LUM Run 3 page) |

AN-2016/325 used the 2016 ReReco value **35.9 fb&#8315;&sup1;**; the UL 2016 total is
36.31 fb&#8315;&sup1;.

---

## 4. Luminosity uncertainty (datacard nuisance)

- **AN-2016/325 (2016)**: `lumi` `lnN` = **2.6%** (`cms-systematics-statistics/references/systematics.md` &sect;2).
- **Run 2 UL**: use the per-year values and the `lumi_13TeV_2016/2017/2018` +
  `lumi_13TeV_correlated` / `lumi_13TeV_1718` correlation scheme from **P1**
  (`LumiRecommendationsRun2`). Not transcribed here — **[Verify]** for the era set
  in the fit.
- A single 2016 fit &rarr; the 2.6% flat value is acceptable as a starting point;
  a 2016+2017+2018 combination needs the P1 correlated / uncorrelated split.

---

## 5. brilcalc (P2) — recomputing the number

To integrate the luminosity actually analysed (e.g. from a CRAB
`processedLumis.json`, or the JSON written by `LumiDumper.py`):

```bash
source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json \
  -u /fb -i <your.json>
```

- `--normtag normtag_PHYSICS.json` is **mandatory** for a physics number — without
  it you get the (unstable) online luminosity.
- `-b "STABLE BEAMS"`; `-o out.csv` for a CSV dump.

---

## 6. Review checklist

1. Era identified; **data** gets the golden JSON (&sect;2), **MC** gets the
   `config` `lumi:` scale (&sect;3).
2. The `data/golden_json/` file `post_proc.py` selects matches the standard UL
   certification for the year.
3. MC luminosity scale matches the LUM UltraLegacy value for the era.
4. Datacard `lumi` nuisance: single 2016 &rarr; 2.6% flat (AN-2016/325); multi-year
   &rarr; the P1 `lumi_13TeV_*` correlation scheme.
5. Any 2022 luminosity value / uncertainty / golden JSON confirmed against the
   current LUM Run 3 page.

---

## 7. Evidence summary

| Item | POG / source | Eras | Established? |
|------|--------------|------|--------------|
| Golden JSON applied to data via `jsonInput` | implementation | all | yes |
| Golden-JSON file names (UL 2016/2017/2018) | DQM-DC | Run 2 UL | yes (standard) |
| MC luminosity scale from `config` `lumi:` | implementation | all | yes; values match LUM UL |
| 2016 luminosity uncertainty 2.6% | HIG / AN-2016/325 | 2016 | yes (2016) |
| Run 2 UL per-year uncertainty + `lumi_13TeV_*` correlation | LUM | 2016&ndash;2018 | **use P1 — [Verify]** |
| 2022 luminosity value / uncertainty / golden JSON | LUM | 2022 | **Authoritative CMS verification required** |

## Last verified

- AN-2016/325 &sect;8 transcription + repo cross-check (`post_proc.py`,
  `config/Input_<year>.yml`): skill update.
- Run 2 UL LUM recommendation (P1): **not transcribed — [Verify]**.
