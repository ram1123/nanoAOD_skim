---
name: cms-systematics-statistics
description: Systematic uncertainties and the statistical analysis (shape-based transverse-mass fit, Higgs Combine datacard, upper limits) for the X/H->ZZ->2l2nu search. Use when building or reviewing datacards, nuisance parameters, K-factors, the interference signal model, shape templates (higgs_combine/), or the systematics table.
---

# CMS Systematics &amp; Statistics — X/H&rarr;ZZ&rarr;2&ell;2&nu;

Companion to `cms-object-guidelines`. That skill covers *object* selections and
their per-object corrections/SFs; **this skill covers the analysis-level
systematic-uncertainty model and the statistical extraction** — the `M_T`
shape-based fit, the Higgs-Combine datacard, the signal + interference model, and
the upper-limit procedure.

**Primary source: CMS AN-2016/325** (2&ell;2&nu; analysis note, 13 TeV, 35.9 fb&#8315;&sup1;,
CMSSW_8_0_25). &sect;6 (optimization + shape fit), &sect;8 + Appendix A (systematics),
&sect;9 (results, interference model Eq. 14), transcribed into the reference files.
It is a **legacy 2016** note — the *method* is authoritative; the *numeric sizes*
(2.6% lumi, 15% Top/WW, &hellip;) are 2016 and must be re-derived per era.

**This is a knowledge/reference skill, not a runnable driver.** The only runnable
surface is the standalone ROOT macros under `higgs_combine/` (`root -l -b -q
higgs_combine/make_shapes.cpp`) — those need `cmsenv` for ROOT, read input
ntuples from EOS paths hard-coded in the macro, and do **not** run full Combine
(`text2workspace.py` / `combine`) which is not set up in this repo. See
`references/statistics.md` &sect;5 for their current state (they are an early
scaffold, not the AN's full model).

## When to use

- writing or reviewing a **Combine datacard** or its shape-input ROOT file;
- deciding which **nuisance parameters** to include, their type (`lnN` / `shape`),
  size, and **correlation** across categories / years / processes;
- the **jet-bin (0 / &ge;1 / VBF) categorization** uncertainty (Stewart&ndash;Tackmann);
- **theory K-factors**: q&#773;qZZ NLO EWK (&rho;-based uncertainty), NNLO QCD, WZ EWK;
- the **signal + interference model** (Eq. 14) and MELA reweighting bookkeeping;
- the **`M_T` shape fit**, blinding, expected/observed limits, CL_s;
- the **data-driven background uncertainties** (non-resonant &alpha;-method 15%;
  Z+jets &gamma;+jets 3-component).

## Required context

1. era(s) in the fit and whether they are combined (&rarr; correlation model);
2. which **categories** enter (`=0j` / `&ge;1j` / `VBF`, `ee` / `&mu;&mu;`, ggF / VBF
   production);
3. the **signal hypothesis** (mass, width &Gamma; = 5/10/100 GeV, interference
   scenario);
4. whether a nuisance is **instrumental** (object SF/scale — cross-check against
   `cms-object-guidelines`) or **theoretical** (xsec/scale/PDF/EWK/UEPS);
5. whether the number's only source is **AN-2016/325** &rarr; mark **[Verify]**
   against the current UL / Run 2+3 2&ell;2&nu; note.

## Reference files — read only what the task needs

- `references/systematics.md` — the full systematic-uncertainty model: the AN &sect;8
  list with 2016 sizes, the Appendix-A **Combine nuisance names** and their
  per-process applicability, instrumental vs theoretical split, jet-bin
  Stewart&ndash;Tackmann, the EWK-K-factor &rho;-method, data-driven-background
  uncertainties, and a cross-check vs this repo (`higgs_combine/make_datacard.cpp`).
- `references/statistics.md` — the statistical extraction: the `M_T` observable and
  binning, the shape-based CL_s limit procedure, the signal-strength +
  interference model (Eq. 14), MELA reweighting, blinding, and the
  `higgs_combine/` macro wiring + its gaps vs the AN.

For the **channel selection / categorization definitions themselves** (not their
uncertainties) read `cms-object-guidelines/references/hzz-2l2nu.md`. For a
**per-object** SF or scale uncertainty (lepton ID/iso/scale, JES/JER, b-veto,
pileup, lumi value) read the matching `cms-object-guidelines` reference — this
skill only says *how it enters the datacard*.

## Evidence rules

Same as `cms-object-guidelines`: for every claimed requirement give the source
(AN section / POG / repo file), the era it applies to, and the last verification
date. Never invent a nuisance size or correlation. If it is only in AN-2016/325,
say `Authoritative CMS verification required` for the current era.

## Reporting categories

- datacard / nuisance defect (wrong type, size, sign, or correlation);
- missing systematic (in AN-2016/325 &sect;8, absent from the datacard);
- analysis-specific inconsistency (repo vs AN method);
- optional improvement;
- authoritative verification required (number only from the 2016 note).
