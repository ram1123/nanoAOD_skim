# MELA Discriminants — Analysis Discriminant (H→ZZ→4l)

Responsible: **analysis-specific matrix-element discriminants — no CMS POG
recommendation applies.** This file documents how MELA is wired into the skim so a
review can check internal consistency, not compliance.

> This replaces the H→μμ `ggh-bdt.md` / `vbf-dnn.md` files, which do not apply to
> this analysis. There is no BDT or DNN here; the 4l channel uses JHUGenMELA
> matrix-element kinematic discriminants. The 2l2q / 2l2nu channels currently
> compute no MELA discriminant.

## Stored sources

| # | Source | Location | Snapshot / verified |
|---|--------|----------|---------------------|
| C1 | JHUGenMELA package (built external dependency) | `external/JHUGenMELA/MELA` — clone of `github.com/JHUGen/JHUGenMELA` **tag `v2.4.2`**, patched with `external/JHUGen_py2to3.patch` | this repo |
| C2 | MELA shared libs loaded at runtime | `modules/H4LCppModule.py` `loadLibraries()` — `libmcfm_710.so`, `libJHUGenMELAMELA.so`, `libjhugenmela.so`, `libcollier.so` from `external/JHUGenMELA/MELA/data/el9_amd64_gcc12/` | this repo |
| C3 | Discriminant construction | `src/H4LTools.cc` ZZ 4l block (~L1000–1085); header `include/H4LTools.h` (`Mela *mela`, L644; `new Mela(13.0, 125.0, TVar::SILENT)`, `setCandidateDecayMode(TVar::CandidateDecay_ZZ)`, L721–722) | this repo |
| C4 | g-constant splines | `external/CoupleConstantsForMELA/gConstant_HZZ2e2mu_{g2,g4,L1,L1Zgs}.root` → `H4LTools::getD{g2,g4,L1,L1Zgs}Constant(massZZ)` | this repo |
| C5 | Output branches | `modules/H4LCppModule.py` `beginFile()` / `analyze()` — `D_CP`, `D_0m`, `D_0hp`, `D_int`, `D_L1`, `D_L1Zg` | this repo |

Classification tags: **[Analysis-specific]**, **[Implementation]**, **[Verify]**.

---

## 1. Setup

- One `Mela` instance per `H4LTools` worker: `new Mela(13.0, 125.0, TVar::SILENT)`
  (√s = 13 TeV, mH = 125 GeV), decay mode `CandidateDecay_ZZ`. **√s is hard-coded
  to 13 TeV** — not adjusted for 2022 (13.6 TeV) inputs. **[Verify]** whether that
  matters for the intended use.
- Four `TSpline` g-constants are read once from `external/CoupleConstantsForMELA/`
  in the constructor (files opened by **relative path** → `post_proc.py` must run
  with CWD = the package directory). All four spline files are the `HZZ2e2mu`
  flavour regardless of the actual 4l final state.
- The MELA `.so`s must be on `LD_LIBRARY_PATH` before `post_proc.py` starts
  (`set_env.sh`, or `eval $(external/JHUGenMELA/MELA/setup.sh env)`).

## 2. Discriminants (`src/H4LTools.cc`, filled only when `ZZSelection_4l()` passes)

Each is built from JHUGen/MCFM probabilities via `mela->computeP(...)` /
`computePM4l(...)` after `setInputEvent` on the selected 4-lepton candidate:

| Branch | Formula (schematically) | g-constant |
|---|---|---|
| `D_0m`  | `me_0plus / (me_0plus + p0minus · getDg4Constant(mZZ)²)` | g4 |
| `D_CP`  | `pg1g4 / (2·√(me_0plus · p0minus))` | — |
| `D_0hp` | `me_0plus / (me_0plus + p0plus · getDg2Constant(mZZ)²)` | g2 |
| `D_int` | `p_ghz1_1_ghz2_1 / (2·√(me_0plus · p0plus))` | — |
| `D_L1`  | `me_0plus / (me_0plus + (p_ghz1prime2_1E4/1e8) · getDL1Constant(mZZ)²)` | L1 |
| `D_L1Zg`| `me_0plus / (me_0plus + (p_ghza1prime2_1E4/1e8) · getDL1ZgsConstant(mZZ)²)` | L1Zgs |

Several `computeP` calls that feed `D_int` / anomalous-coupling terms carry
`//FIXME` in the source — treat those discriminants as provisional until checked.

`D_bkg_kin` / `D_bkg` are declared in the header but **not** filled into output
branches by the current module.

## 3. What a review should check

- √s (13.0) and mH (125.0) hard-coded values vs the samples being processed.
- The `HZZ2e2mu` g-constant splines applied to 4e / 4μ candidates — is that the
  intended (mass-only) treatment?
- Relative-path `TFile::Open` for the g-constants (CWD dependency; no existence check).
- The `//FIXME` `computeP` probability configurations for `D_int`, `D_L1*`.
- MELA `setInputEvent` / `resetInputEvent` pairing per event (leak / stale-candidate check).
- That discriminant branches are `-999` sentinels for non-4l events and consumers
  downstream (`higgs_combine/`) filter on `passZZ4lSelection`.
- MELA/JHUGen version pin: `external/JHUGenMELA` at `v2.4.2` — record it in any report.

## 4. Not applicable here

No CMS-POG recommendation governs these discriminants. Do not classify a MELA
finding as an "official recommendation violation" — use **[Analysis-specific]** or
**[Implementation]**.

## Last verified

- Local source review: (fill in on first use)
- Upstream JHUGenMELA: tag `v2.4.2` (as pinned in `docs/README.md` / `setup.sh`)
