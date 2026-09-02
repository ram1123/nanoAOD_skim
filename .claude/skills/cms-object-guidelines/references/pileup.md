# Pileup Reweighting — 2&ell;2&nu; &amp; LUM Recommendations

> For the **X/H&rarr;ZZ&rarr;2&ell;2&nu; NanoAOD-tools skim**
> (`post_proc.py` &rarr; `modules/H4LCppModule.py` &rarr; C++ `src/H4LTools.cc`).
>
> - Eras: **Run 2 UL 2016 / 2017 / 2018** and some **2022**; NanoAOD **v9** and
>   **v15**. Run-3-2023+ material does not apply.
> - PU reweighting: nanoAOD-tools **`puAutoWeight_2016 / 2017 / 2018`** producer,
>   added for **MC only** in `post_proc.py` (the "PU-weight module" in the chain).
>   **None wired for 2022** — **[Verify]** / resolve before running 2022 MC.
> - Legacy Run 2 PU-jet-ID scale factors: `modules/JetSFMaker.py`, only under
>   `post_proc.py --WithSyst` (see `references/jets.md` &sect;6).
> - Re-check the minimum-bias cross section against the LUM recommendation.

Responsible POG: **LUM** (pileup reweighting, minimum-bias cross section).
Channel strategy: `references/hzz-2l2nu.md`.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| C1 | This repo | nanoAOD-tools `puAutoWeight_2016/2017/2018` producer in `post_proc.py` (MC only); `modules/JetSFMaker.py` (PU-jet-ID SF, `--WithSyst`) | repo |
| P1 | LUM entry point | `twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData` / `LumiRecommendationsRun2` (`#Pileup_reweighting`); LUM `puWeights.json.gz` (`jsonpog-integration` `POG/LUM`) | reference |

Not covered by stored sources &rarr; **Authoritative CMS verification required**:
the minimum-bias cross section and its &plusmn; variation for the target era, the
2022 pileup-reweighting prescription/payload, and the LUM `puWeights.json.gz`
version.

Classification tags: **[LUM official]**, **[Implementation]** (this repo),
**[Verify]**.

---

## 1. Pileup reweighting (C1)  **[LUM official]**

- An event weight correcting the MC in-time pileup profile (`Pileup_nTrueInt`) to
  the data profile. Applied to **simulation only**.
- This repo uses the nanoAOD-tools **`puAutoWeight_20XX`** module for Run 2 UL
  (2016 / 2017 / 2018), which pulls the LUM data pileup profile and computes the
  weight internally.
- **2022**: no `puAutoWeight` is wired in `post_proc.py` — **[Verify]** which LUM
  Run 3 `puWeights.json.gz` / minimum-bias value to use and add the producer
  before a 2022 MC run.

### 1.1 Minimum-bias cross section &amp; variation

- Nominal `&sigma;_minbias = 69.2 mb`; the up / down pileup weights come from
  shifting it by `&plusmn; ~4.6%` (i.e. `~66.0` / `~72.4 mb`).
- **[Verify]** the exact value and variation against the LUM recommendation for
  the era — do not assume.
- Carried as a **weight systematic** (`pileup` up / down), one of the
  `cms-systematics-statistics` shape nuisances (`references/systematics.md` &sect;2).

---

## 2. Legacy Run 2 PU jet ID

For Run 2 **CHS** (v9) jets, PU jet ID matters for `pT < 50 GeV`; its scale
factors are applied only via `modules/JetSFMaker.py` under `--WithSyst`. Working
points and SFs: the JME PU jet ID TWiki (`references/jets.md` &sect;6). Run 2 /
Run 3 **PUPPI** (v15) jets do not use a legacy PU jet ID — pileup is handled by
the PUPPI weights.

---

## 3. Review checklist

1. PU reweighting applied to **MC only**; `puAutoWeight_20XX` present for the Run 2
   UL year; 2022 prescription resolved.
2. Minimum-bias cross section (69.2 mb) and its &plusmn; variation confirmed for the
   era; `pileup` up / down wired into the weight systematics.
3. Run 2 v9: if `--WithSyst`, PU-jet-ID SF applied via `JetSFMaker`
   (`references/jets.md` &sect;6). Run 2/3 v15 PUPPI: no legacy PU jet ID.

---

## 4. Evidence summary

| Item | POG / source | Eras | Established? |
|------|--------------|------|--------------|
| PU reweighting on MC via `puAutoWeight_20XX` | LUM / implementation | 2016&ndash;2018 | yes |
| PU reweighting for 2022 | LUM | 2022 | **not wired** — **[Verify]** |
| Minimum-bias xsec (69.2 mb) + variation | LUM | all | value standard; **[Verify]** per era |
| `pileup` up/down weight systematic | LUM | all | yes (`cms-systematics-statistics`) |
| Legacy Run 2 PU-jet-ID SF (`JetSFMaker`, `--WithSyst`) | JME | 2016&ndash;2018 v9 | yes (implementation) |

## Last verified

- Repo cross-check (`post_proc.py`, `modules/JetSFMaker.py`): skill update.
- Current LUM minimum-bias / Run 3 pileup recommendation: **not consulted — [Verify]**.
