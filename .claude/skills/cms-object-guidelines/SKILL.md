---
name: cms-object-guidelines
description: Apply official CMS POG recommendations for object selection, corrections, scale factors, and uncertainties in the X/H->ZZ->2l2nu NanoAOD skim.
---

# CMS Object Guidelines — X/H&rarr;ZZ&rarr;2&ell;2&nu;

Use this skill for work involving:

- the **2l2nu channel strategy** — selection chain, transverse mass, event
  categorization (0-jet / &ge;1-jet / VBF), data-driven backgrounds, K-factors,
  and the systematic-uncertainty list;
- muons;
- electrons;
- FSR-photon recovery for leptons;
- AK4 jets (categorization, VBF tag jets, b-veto);
- b tagging (the 2l2nu b-veto and the &alpha;-method b-tag requirement);
- missing transverse momentum and MET filters — **2l2nu is MET-driven**: the
  `MET > 125 GeV` cut and the d&Phi;(jet,MET) / d&Phi;(Z,MET) cuts define the
  signal region;
- MET-&phi; (XY) correction;
- trigger selection and trigger-object matching;
- object corrections and scale factors;
- object-related systematic uncertainties;
- integrated luminosity and golden-JSON lumimask;
- pileup reweighting.

The MELA matrix-element discriminants (`D_CP`, `D_0m`, …) belong to the H&rarr;ZZ&rarr;4&ell;
path, not 2&ell;2&nu; (which fits `M_T` directly); `references/mela-discriminants.md`
documents them only for internal-consistency checks.

## Analysis context — read first

**Focus: X/H&rarr;ZZ&rarr;2&ell;2&nu;** (branch `HZZ_Analysis_2l2nu_dev`, NanoAOD v15;
also v9). When a request does not name a channel, assume **2l2nu**. The
channel-strategy reference is `references/hzz-2l2nu.md`, transcribed from CMS
**AN-2016/325** (the 2&ell;2&nu; analysis note — 2016/legacy: the *method* is
authoritative, the *numeric working points* are superseded and marked **[Verify]**).

**Sibling skill:** `cms-systematics-statistics` covers the analysis-level
systematic-uncertainty model and the statistical extraction (the `M_T` shape fit,
the Higgs Combine datacard, nuisance names, K-factor uncertainties, the Eq. 14
interference signal model, `higgs_combine/`). This skill stays at the **object**
level — it says what a correction / SF *is*; the sibling says how it *enters the
datacard*.

### The skim

- `post_proc.py` builds a `nanoAOD-tools` `PostProcessor` chain; the physics
  module is `modules/H4LCppModule.py` (`HZZAnalysisCppProducer`), which drives the
  C++ worker `src/H4LTools.cc` / `include/H4LTools.h`.
- **All cut values come from `config/Input_<year>.yml`**, pushed into the worker
  via `H4LTools::Initialize*` setters. Change thresholds / triggers there; change
  selection logic in `src/H4LTools.cc`.
- Eras in scope: **Run 2 UL 2016 / 2017 / 2018** (13 TeV) and some **2022**
  (13.6 TeV). NanoAOD **v9** and **v15** (`config/Input_2016.yml` is referenced
  but absent). Run-3-2023+ material does not apply.
- Muons and electrons are both **primary** selection objects (they build the
  Z&rarr;&mu;&mu; / Z&rarr;ee leg); **FSR recovery** is applied to both; MET is the
  Z&rarr;&nu;&nu; leg. AK4 jets give the jet-multiplicity category, the VBF tag
  jets, and the b-veto.

The `references/*.md` files carry the official POG recommendation for each object
as a **starting point**; re-verify every number against the current 2&ell;2&nu;
analysis note / HIG working group before treating it as a requirement. A value
whose only source is AN-2016/325 (2016/legacy) is **[Verify]** for the target era.

## Required context

Before applying a recommendation, identify:

1. era (2016preVFP / 2016postVFP / 2017 / 2018 / 2022);
2. data or simulation;
3. NanoAOD version (v9 &rarr; `Electron_mvaFall17V2Iso_WP90`, stored `Jet_jetId`;
   v15 &rarr; `Electron_mvaIso_WP90`, jet ID recomputed by `H4LTools::PassJetIDv15`);
4. path — **default 2l2nu signal region**; also the 2l2nu e&mu; control region;
5. intended working point;
6. responsible CMS POG (and, for a channel-strategy question, `references/hzz-2l2nu.md`).

If something material is unknown, ask or mark the conclusion unverified.

## Reference selection

Read only the relevant file(s):

- `references/hzz-2l2nu.md` — **2l2nu channel strategy** (AN-2016/325):
  signal-region selection, e&mu; control region, transverse mass, categorization,
  background methods, K-factors, systematics; includes a cross-check table vs the
  current skim. Read this first for any 2l2nu selection / categorization /
  background question, then the object file for the specific object.
- `references/muons.md` — Z&rarr;&mu;&mu; leg, Rochester (`muonScaleRes20XX`), FSR
- `references/electrons.md` — Z&rarr;ee leg, EGM MVA-iso WP90 / HZZ BDT points, FSR
- `references/jets.md` — AK4 jets: JEC/JER (`--WithSyst`), v15 jet-ID recompute,
  jet veto maps, the 2l2nu jet cuts / VBF tag
- `references/b-tagging.md` — DeepJet b-veto
- `references/met.md` — PuppiMET, MET-&phi; correction, MET filters, the 2l2nu
  MET-driven selection
- `references/lumi.md` — integrated luminosity, golden JSON, lumi uncertainty
- `references/pileup.md` — pileup reweighting (`puAutoWeight_20XX`)
- `references/mela-discriminants.md` — JHUGenMELA discriminants (4&ell; only;
  internal-consistency checks, not a 2l2nu requirement)

Do not load every reference automatically. Trigger-object matching and
overlap/cleaning are covered inside the object files. Photons (as prompt objects)
and taus are not used — treat such a request as `Authoritative CMS verification
required`.

## Review checklist

When inspecting an object implementation, verify:

- kinematic acceptance (`config/Input_<year>.yml` vs `src/H4LTools.cc`);
- identification working point (and the v9 / v15 branch name);
- isolation definition and working point;
- impact-parameter requirements;
- FSR-recovery inputs and selection;
- data-quality / MET filters (`modules/METFilters.py`);
- correction sequence (Rochester `muonScaleRes20XX`, MET-&phi;, JEC/JER under
  `--WithSyst`);
- data/MC and trigger scale factors (the nominal skim applies **none** for
  leptons / b-tag);
- object cleaning and overlap-removal order;
- systematic-uncertainty variations;
- applicability to the era and NanoAOD version in use;
- cut-flow bin labels (`dynamicCuts_*` in `H4LCppModule`) matching the mask
  actually applied.

For **2l2nu** work also check, against `references/hzz-2l2nu.md` &sect;2 / &sect;9:

- the dilepton mass window `|m_ll − 91| < 15` and `p_T^ll > 55` are actually
  applied (config currently has `M_ll_Window: 0`, `Pt_ll: 10`);
- muon `|η| < 2.4` in the 2l2nu block (config uses `Lep_eta: 2.5` for both flavours);
- the b-jet veto is active (currently commented out in `ZZSelection_2l2nu()`);
- the final `MET > 125 GeV` cut rejects events (code only counts `> 100`);
- `|Δφ(Z, MET)| > 0.5` is present;
- the transverse-mass definition matches note Eq. 6 (Z→νν leg at `m_Z`);
- VBF centrality + central-jet veto, and the qqZZ EWK / NNLO-QCD K-factors.

## Evidence rules

For every claimed official requirement, provide:

- responsible POG;
- applicable era;
- source URL or official repository;
- stored source version or Git tag, when available;
- last verification date.

Never invent or extrapolate an official recommendation. If the stored
documentation is incomplete, or a number's only source is AN-2016/325
(2016/legacy) and it has not been re-verified for the target era, say:

`Authoritative CMS verification required.`

## Reporting categories

Classify findings as:

- official recommendation violation;
- analysis-specific inconsistency;
- implementation defect;
- optional improvement;
- authoritative verification required.
