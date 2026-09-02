---
name: cms-object-guidelines
description: Apply stored official CMS POG recommendations for object selection, corrections, scale factors, and uncertainties.
---

# CMS Object Guidelines

Use this skill for work involving:

- the **2l2nu channel strategy** — selection chain, transverse mass, event
  categorization (0-jet / ≥1-jet / VBF), data-driven backgrounds, K-factors, and
  the systematic-uncertainty list;
- muons;
- electrons;
- FSR-photon recovery for leptons;
- AK4 or AK8 jets (resolved and boosted Z→qq tagging);
- b tagging (2l2nu b-veto, 2l2q);
- missing transverse momentum and MET filters (**2l2nu is MET-driven** — the
  MET > 125 GeV cut and dΦ(jet,MET) / dΦ(Z,MET) cuts define the signal region);
- MET-φ (XY) correction;
- trigger selection and trigger-object matching;
- object corrections and scale factors;
- object-related systematic uncertainties;
- integrated luminosity and golden-JSON lumimask;
- pileup reweighting;
- the MELA kinematic discriminants (`D_CP`, `D_0m`, …) used in the 4l channel.

## Analysis context — read first

**Current focus: X/H→ZZ→2l2nu.** The active work targets the **2l2nu** channel
only (branch `HZZ_Analysis_2l2nu_dev`, NanoAOD v15). The 4l and 2l2q code paths
still exist but are secondary — parts are commented out on this branch. When a
request does not name a channel, assume **2l2nu**. The channel-strategy reference
is `references/hzz-2l2nu.md` (transcribed from CMS **AN-2016/325**, the 2l2nu
analysis note — 2016/legacy, method is authoritative, numeric WPs superseded).

This repository is the **H→ZZ→4l / 2l2q / 2l2nu NanoAOD skim**:

- `post_proc.py` builds a `nanoAOD-tools` `PostProcessor` chain; the physics module
  is `modules/H4LCppModule.py` (`HZZAnalysisCppProducer`), which drives the C++
  worker `src/H4LTools.cc` / `include/H4LTools.h`.
- **All cut values come from `config/Input_<year>.yml`**, pushed into the worker
  via `H4LTools::Initialize*` setters. Change thresholds/triggers there; change
  selection logic in `src/H4LTools.cc`.
- Eras in scope: **Run 2 UL 2016 / 2017 / 2018** and some **2022**. NanoAOD **v9**
  and **v15** (`config/Input_2016.yml` is referenced but absent).
- Object roles differ from a single-final-state analysis: **muons and electrons
  are both primary selection objects** (4e / 4μ / 2e2μ, and Z→ll in 2l2q/2l2nu);
  **FSR recovery** is applied to both; **2l2nu is MET-driven**; AK8 fat jets tag
  the boosted Z→qq in 2l2q.

The `references/*.md` files were **ported from the H→μμ (copperhead) analysis** and
carry a banner to that effect. Their CMS-POG recommendation content is a starting
point only — the coffea `src/copperhead_processor.py` / `configs/parameters/*.yaml`
/ stage-1/stage-2 machinery they cite **does not exist here**, and Run-3-2023+
material does not apply. Re-verify against this analysis's H→ZZ note / HIG group
before treating anything as a requirement.

## Required context

Before applying a recommendation, identify:

1. era (2016 / 2016APV / 2017 / 2018 / 2022);
2. data or simulation;
3. NanoAOD version (v9 → `Electron_mvaFall17V2Iso_WP90`, stored jet ID/PU-ID;
   v15 → `Electron_mvaIso_WP90`, jet ID recomputed by `H4LTools::PassJetIDv15`);
4. channel (**default 2l2nu**; also 2l2nu-emu-CR / 4l / 2l2q) — it changes which
   objects matter and which reference to open;
5. intended working point;
6. responsible CMS POG (and, for a channel-strategy question, `references/hzz-2l2nu.md`).

If something material is unknown, ask or mark the conclusion unverified. A number
whose only source is AN-2016/325 (2016/legacy) is **[Verify]** until confirmed
against the current UL / Run 2+3 2l2nu note.

## Reference selection

Read only the relevant file(s):

- `references/hzz-2l2nu.md` — **2l2nu channel strategy** (AN-2016/325): signal-region
  selection, e_mu control region, transverse mass, categorization, background
  methods, K-factors, systematics; includes a cross-check table vs the current
  skim. Read this first for any 2l2nu selection / categorization / background
  question, then the object file for the specific object.
- `references/muons.md`
- `references/electrons.md`
- `references/jets.md` — AK4 jets: JEC/JER, jet ID (incl. v15 recompute), PU jet ID
- `references/fat-jets.md` — AK8 jets and substructure (boosted 2l2q)
- `references/b-tagging.md`
- `references/met.md` — includes MET noise / event filters and MET-φ correction
- `references/lumi.md` — integrated luminosity, golden JSON, lumi uncertainty
- `references/pileup.md` — pileup reweighting
- `references/mela-discriminants.md` — JHUGenMELA discriminants (analysis-specific,
  no CMS-POG recommendation; use it for internal-consistency checks only)

Do not load every reference automatically. Trigger-object matching and
overlap/cleaning are covered inside the object files. Photons (as prompt objects)
and taus are not used — treat such a request as `Authoritative CMS verification
required`.

## Review checklist

When inspecting an object implementation, verify:

- kinematic acceptance (`config/Input_<year>.yml` vs `H4LTools`);
- identification working point (and the v9/v15 branch name);
- isolation definition and working point;
- impact-parameter requirements;
- FSR-recovery inputs and selection;
- data-quality / MET filters;
- correction sequence (Rochester `muonScaleRes*`, MET-φ, JEC/JER under `--WithSyst`);
- scale and resolution corrections; data/MC and trigger scale factors;
- object cleaning and overlap-removal order;
- systematic uncertainty variations;
- compatibility between selection and scale factors;
- applicability to the era and NanoAOD version in use;
- cut-flow bin labels (`dynamicCuts_*` in `H4LCppModule`) matching the mask actually applied.

For **2l2nu** work also check, against `references/hzz-2l2nu.md` §2 / §9:

- the dilepton mass window `|m_ll − 91| < 15` and `p_T^ll > 55` are actually
  applied (config currently has `M_ll_Window: 0`, `Pt_ll: 10`);
- muon `|η| < 2.4` in the 2l2nu block (config uses `Lep_eta: 2.5` for both flavours);
- the b-jet veto is active (currently commented out in `ZZSelection_2l2nu()`);
- the final `MET > 125 GeV` cut rejects events (code only counts `> 100`);
- `|Δφ(Z, MET)| > 0.5` is present;
- the transverse-mass definition matches note Eq. 6 (Z→νν leg at `m_Z`);
- VBF centrality + central-jet veto, and the qqZZ EWK/NNLO-QCD K-factors.

## Evidence rules

For every claimed official requirement, provide:

- responsible POG;
- applicable era;
- source URL or official repository;
- stored source version or Git tag, when available;
- last verification date.

Never invent or extrapolate an official recommendation. If the stored
documentation is incomplete or was ported from H→μμ without re-verification, say:

`Authoritative CMS verification required.`

## Reporting categories

Classify findings as:

- official recommendation violation;
- analysis-specific inconsistency;
- implementation defect;
- optional improvement;
- authoritative verification required.
