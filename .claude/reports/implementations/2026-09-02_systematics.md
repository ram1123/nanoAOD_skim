# Systematic-uncertainty review & weight-branch implementation

- Date: 2026-09-02
- Type: Implementation
- Status: Partially resolved (weight-only systematics done; shape systematics + SF weights deferred)
- Applicable era: Run 2 UL 2016preVFP / 2016postVFP / 2017 / 2018
- NanoAOD campaign: v9 and v15
- Relevant files: `post_proc.py`, `modules/H4LCppModule.py`, `modules/keep_and_drop_list.py`, `docs/systematics_review.md`

## Question

Full systematic-uncertainty review of the X&rarr;ZZ&rarr;2l2nu skim: what is
implemented, is nominal&harr;shifted consistent, what CMS-recommended systematics
are missing, implement the clearly-required ones, add nothing irrelevant. Run
through `physics-reviewer` + `code-reviewer` + `test-specialist` via `/coordinate`.

## Evidence

- Reviewed `docs/systematics_review.md` (an earlier main-session pass) against the
  code as it stood.
- `physics-reviewer`: confirmed most rows; found the pileup producer targets the
  pre-UL data profile, the "2018 prefire &rarr; nominal" claim is physically
  wrong (muon prefiring is real for 2018), lepton/b-tag SFs are pure weights not
  re-selection items, and qqZZ/ggZZ K-factors are nominal corrections that are
  currently un-applicable (gen kinematics dropped).
- `code-reviewer`: all 11 plumbing checks pass; independently flagged the pileup
  producer campaign, the downstream double-count risk (`make_shapes.cpp` uses
  bare `puWeight` + `(nMC-2nNeg)`), and the hard-coded `AK4PFchs` in `--WithSyst`.
- Confirmed in `treeReaderArrayTools.py:80` that an unknown branch raises
  `RuntimeError`, which `getattr(event, name, default)` does **not** catch &mdash;
  so the Pass-1 weight code would crash (not silently fall back) for a year with
  no PU / prefiring producer (2022).
- `test-specialist` smoke test (2l2nu, 2018 v15 signal, `-n 300`, inside
  `cmssw-el9`): exit 0, 166 selected. All nine weight branches present.
  `overallEventWeight` varies 166/166 and matches
  `genWeight&middot;puWeight&middot;L1PreFiringWeight_Nom` on entry 0. PU up&ne;down
  166/166; prefire up&ne;down 92/166 (rest have `L1PreFiringWeight_*=1`). No
  fallback warnings, no NaN/Inf, no traceback.

## Findings

Confirmed:
- `overallEventWeight` plumbing (ordering, reachability, branch parity, keep
  rules, no code-level double application, CRLF) is correct.
- PU and L1-prefiring are the only weight systematics that need no re-selection
  and were the correct clearly-required subset to implement.

Corrected during review:
- Pileup producer must be `puWeight_UL20XX` (UL data profile + fixed MC profile),
  not `puAutoWeight_20XX` (pre-UL ReReco profile, per-file MC profile).
- `getattr(event, name, default)` is not a valid guard against a missing branch.
- 2018 L1-prefiring (muon term) is non-negligible; up/down do not collapse.

Assumptions / open:
- `puWeight_UL2016` has no preVFP/postVFP split; 2022 has no PU payload. Handled
  by a one-time warning + `puWeight` left out of `overallEventWeight`.
- 2016 is non-functional on this branch (`config/Input_2016.yml` absent, no v15
  token, no APV split).

## Decision or outcome

Implemented (weight-only, no selection change):
- `post_proc.py`: `puAutoWeight_{2016,2017,2018}` &rarr; `puWeight_UL{...}`;
  `--WithSyst` jetType `AK4PFPuppi` for v15 else `AK4PFchs` (MC + data blocks).
- `modules/H4LCppModule.py`: `_get_weight_branch()` helper (try/except, one-time
  warning); `analyze()` weight block uses it; `beginFile()` comment expanded into
  a downstream normalisation contract (normalise by `genEventSumw`; do not
  re-apply `genWeight`/`puWeight`; fix `make_shapes.cpp` accordingly).
- Branches (from Pass 1): `overallEventWeight`, `overallEventWeight_pu{Up,Down}`,
  `overallEventWeight_prefire{Up,Down}` &mdash; MC only, data = 1.0.
- `docs/systematics_review.md` rewritten with both reviews folded in.

Deferred (needs external payload &mdash; implementable as a weight, no re-selection):
muon/electron reco/ID/iso SF, dilepton trigger SF, b-tag SF (also needs the
b-veto re-enabled + a DeepJet payload; bundled CSVs are DeepCSV),
qqZZ/ggZZ EWK+QCD K-factors (also needs gen ZZ kinematics kept).

Deferred (needs a per-variation re-run of the C++ selection): JES, JER,
unclustered-MET, electron energy scale/smearing, Rochester propagation.

Downstream / datacard only: luminosity `lnN`, MC stats (`autoMCStats`),
non-resonant & Z+jets data-driven, LHE scale/PDF/PS envelopes.

## Verification

- Command: `python3 -m py_compile post_proc.py modules/H4LCppModule.py modules/keep_and_drop_list.py`
- Result: clean.
- Command: `python3 post_proc.py -i config/ExampleInputFileList.txt -n 300` (inside `cmssw-el9`)
- Result: exit 0; 166 selected; 9 weight branches present; `overallEventWeight`
  varies and matches the definition; PU & prefire up&ne;down; no warnings/NaN/traceback.
- CRLF of `post_proc.py` preserved (226/226).

## CMS sources

- LUM: UltraLegacy pileup, 69.2 mb min-bias, &plusmn;4.6%; `puWeights.json.gz`
  (`Collisions1{6,7,8}_UltraLegacy_goldenJSON`).
- L1 prefiring: ECAL (2016&ndash;2017) + muon (2016&ndash;2018), `L1PreFiringWeight_*`.
- JME: regrouped JES (~11 sources), JER; PUPPI-MET &phi; correction applicability
  is `[Verify]` for UL.
- BTV: DeepJet `btagging.json.gz`, shape/`comb` method.
- AN-2016/325 &sect;8 + Appendix A (nuisance names) &mdash; 2016 pre-legacy,
  numbers unverified for UL.

## Remaining work

- Add lepton (reco/ID/iso/trigger) and b-tag SF weights once POG JSONs + verified
  WP keys are in hand.
- Update `higgs_combine/make_shapes*.cpp` to the `overallEventWeight` +
  `genEventSumw` contract (drop bare `puWeight`, change the denominator).
- Architect a per-variation selection path in `H4LTools` for JES/JER/unclustered
  MET / lepton scale, with PuppiMET shifts.
- Keep gen ZZ/Z kinematics (or apply the K-factor in `H4LCppModule`) so the
  qqZZ/ggZZ EWK+QCD K-factors become applicable.
- Create `config/Input_2016.yml` + a 2016 preVFP/postVFP split (year detect + PU).
- 2022: LUM Run-3 `puWeights.json.gz`.
