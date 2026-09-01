# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A NanoAOD skimmer for H&rarr;ZZ&rarr;4l / 2l2q / 2l2nu analyses. It is a CMSSW
package built on the `nanoAOD-tools` `PostProcessor` framework and must be checked
out at this exact path:

```
$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
```

Almost everything resolves paths relative to this directory or `$CMSSW_BASE`, so
run commands from here after `cmsenv`.

**This branch (`HZZ_Analysis_2l2nu_dev`) is the 2l2nu working branch and targets
NanoAOD v15.** Compared with `HZZ_Analysis` it defaults to the 2l2nu channel, adds
MET-phi correction, migrates object branch names to v15, and adds standalone
Higgs-Combine macros. The 4l and 2l2q code paths still exist but parts are
commented out here (e.g. FatJet handling in `modules/H4LCppModule.py`).

### Branch / remote layout in this checkout

- `HZZ_Analysis` &mdash; leave untouched unless asked; different (older) analysis state.
- `HZZ_Analysis_2l2nu` &mdash; tracks `origin/HZZ_Analysis_2l2nu` (ram1123 repo).
- `HZZ_Analysis_2l2nu_dev` &mdash; this branch; `HZZ_Analysis_2l2nu` + merged work
  from the `anusree` remote (`github.com/Anusreevijay769/nanoAOD_skim`).
- `origin` uses an SSH URL (`git@github.com:ram1123/nanoAOD_skim.git`). If no SSH
  key is available, fetch over HTTPS, e.g.
  `git fetch https://github.com/ram1123/nanoAOD_skim.git <branch>:refs/remotes/origin/<branch>`.

## Build / setup

Full instructions: [docs/README.md](docs/README.md); [setup.sh](setup.sh) scripts
it end to end. The parts that matter:

- CMSSW release `CMSSW_14_0_2` (`el9_amd64_gcc12`).
- `PhysicsTools/NanoAODTools` is a **fork** (`ram1123/nanoAOD-tools`, branch
  `h4l_allCh_dev`), not upstream. This branch depends on that fork providing
  `postprocessing/modules/common/met_phi_correction.py` (`METPhiCorrector`,
  `Campaign`).
- `external/yaml-cpp` &mdash; separate clone, patched with
  `external/yamlcpp_pkg_py2to3.patch`, built as a shared lib into
  `external/yaml-cpp/build/` (`cmake3 .. -DBUILD_SHARED_LIBS=ON && cmake3 --build .`).
- `external/JHUGenMELA` &mdash; separate clone at tag `v2.4.2`, patched with
  `external/JHUGen_py2to3.patch`, built with `JHUGenMELA/MELA/setup.sh -j 8`;
  `chmod +x` the `.so`s under `MELA/data/el9_amd64_gcc12/`.
- b-tag CSVs: `cp data/btag/*.csv $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF/`.
- `scram b -j 8` from `$CMSSW_BASE/src`.

`external/JHUGenMELA/` and `external/yaml-cpp/` are git-ignored (own repos).
`src/H4LTools.cc` is **not** compiled by `scram`; ROOT ACLiC JIT-compiles it at
runtime (see Architecture).

## Running interactively

```bash
cmsenv
source set_env.sh          # MELA env (LD_LIBRARY_PATH + `eval $(.../MELA/setup.sh env)`) + voms-proxy-init
python3 post_proc.py -i <file.root | list.txt> -n 1000
```

`set_env.sh` runs `voms-proxy-init` and `exit 1`s if the proxy file is missing. If
you only need the environment, run the two MELA lines from `set_env.sh` yourself.

There is **no test suite**. To smoke-test a change, run `post_proc.py` over one
file with a small `--entriesToRun` and `--DEBUG`.

### `post_proc.py` arguments (note the defaults on this branch)

- `-i/--inputFile` &mdash; a `.root` file, or a `.txt` with one LFN per line (each
  gets `root://cms-xrd-global.cern.ch/` prepended). Empty &rarr; reads
  `ExampleInputFileList.txt`.
- `-n/--entriesToRun` &mdash; **default `0` = all events** (was 100 on `HZZ_Analysis`).
- `--channels {all,4l,2l2q,2l2v}` &mdash; **default `2l2v`**.
- `--WithSyst` &mdash; enables JME correctors + PU-ID SF. **Currently broken for MC**:
  `post_proc.py` references undefined `muonScaleRes()` / `gammaSF()` in the
  `if isMC and args.WithSyst` block &mdash; it raises `NameError`. Leave `--WithSyst`
  off, or fix that line, before running systematics.
- `-o/--outputFile`, `-outDir/--outputDir`, `-c/--cutFlowFile`, `--DEBUG`.

**Year / data detection is substring-matching on the first input file's path**
(`post_proc.py` `main()`), and on this branch the tokens are NanoAOD-version
specific and the checks are independent `if`s (not `elif`), so a later match
overrides an earlier one:

| token(s) in path | year / config |
| --- | --- |
| `Summer22`, `Run2022` | 2022 &rarr; `config/Input_2022.yml` |
| `UL18NanoAODv9`, `UL2018_MiniAODv2_NanoAODv9` | 2018 &rarr; `config/Input_2018.yml` |
| `UL17NanoAODv15`, `UL2017` | 2017 &rarr; `config/Input_2017.yml` |
| `20UL16NanoAODAPVv9` / `20UL16NanoAODv9` | 2016 &rarr; `config/Input_2016.yml` (**missing from repo**) |
| `UL2018_NanoAODv15`, `UL18NanoAODv15` | 2018 &rarr; `config/Input_2018.yml` |

`isMC = "/data/" not in first_file`. If no token matches, `year`/`cfgFile` stay
`None` and the job fails.

## Batch submission (HTCondor, lxplus)

```bash
python3 scripts/condor/condor_setup_lxplus.py \
        --input_file input_data_files/sample_list_v15_2018.dat \
        --submission_name HZZ2l2nu_<date> --condor_queue tomorrow \
        --condor_file_name submit_condor_HZZ2l2nu_<date>
voms-proxy-init -voms cms --valid 200:00
condor_submit <generated>.jdl
```

`condor_setup_lxplus.py` reads DAS dataset names (one per line, `#` = comment)
from the path given to `--input_file` **as-is relative to CWD** (on this branch it
is no longer auto-prefixed with `input_data_files/` &mdash; pass the full relative
path). It expands datasets with `dasgoclient`, tars the whole CMSSW area and
`xrdcp`s it to EOS, and writes a `.jdl` + `.sh` + `.txt` triple. It imports
helpers (`color_style`, `infoCreaterGit`, `fileshelper`, `makeTarFile`) from
`scripts/utils/` &mdash; that dir must be importable. Resubmit failures with
`scripts/condor/nanoAOD_condor_resubmit.py`. `scripts/crab/` is an older,
not-recently-tested path.

## Architecture

### Flow

`post_proc.py` assembles an ordered list of `nanoAOD-tools` `Module`s and hands
them to `PostProcessor`: the year's muon scale/resolution producer, the main
analysis module `H4LCppModule`, and (MC only) `GenVarsProducer` + the PU-weight
module; `--WithSyst` would add the JME correctors + `JetSFMaker`. Output branch
selection is built at runtime from the Python lists in
`modules/keep_and_drop_list.py`, written to a temp keep/drop file.
`fwkJobReport=True` is set deliberately so `haddnano.py` runs.

### The analysis module &mdash; two layers

`modules/H4LCppModule.py` (`HZZAnalysisCppProducer`) is a thin Python `Module`:

1. `loadLibraries()` &mdash; loads the MELA `.so`s and `libyaml-cpp.so`, then
   `ROOT.gROOT.ProcessLine(".L <base>/src/H4LTools.cc+O")` to JIT-compile the C++
   worker.
2. `__init__` &mdash; parses `config/Input_<year>.yml` and pushes every cut value
   into the C++ worker via `Initialize*` setters (`InitializeElecut`,
   `InitializeMucut`, `InitializeJetcut`, `InitializeEvtCut`,
   `InitializeHZZ2l2qCut`, `InitializeHZZ2l2nuCut`). Builds the cut-flow `TH1F`
   whose bins come from the `dynamicCuts_*` name lists (the `dynamicCuts_2l2nu`
   list on this branch includes fine-grained per-muon bins like `cut_mu_pt`).
3. `analyze(event)` per event &mdash; resets the worker; copies
   `Electron`/`Muon`/`Jet`/`FsrPhoton`/`MET`/`PuppiMET`/`GenPart`/`GenJet`
   collections into it with `Set*` calls; runs the MET-phi correction (see below);
   calls selection entry points &mdash; `GetZ1_2l2qOR2l2nu()`, `GetZ1_emuCR()`,
   `ZZSelection_4l()`, `ZZSelection_2l2q()`, `ZZSelection_2l2nu()` &mdash; and reads
   results back off worker attributes (`Z1`, `Z2`, `ZZsystem`, `HZZ2l2nu_*`,
   per-cut counters, MELA `D_*`, &hellip;) to fill output branches.

`src/H4LTools.cc` + `include/H4LTools.h` hold the physics: object ID, FSR
recovery, ZZ candidate construction for all channels + the e&mu control region,
and MELA discriminants. Cut thresholds are **not** hard-coded here &mdash; they
arrive from the YAML via `Initialize*`. Change selection logic here; change
thresholds/triggers in the YAML.

### NanoAOD v15 migration (this branch)

`H4LCppModule.analyze` and `H4LTools` were changed for v15 inputs; the v9 code is
left commented alongside:

- Electrons: `SetElectrons(pt, eta, phi, mass, dxy, dz, pdgId, mvaIso_WP90,
  pfRelIso03_all)` &mdash; v15 uses `Electron_mvaIso_WP90` (v9:
  `mvaFall17V2Iso_WP90`), and the argument order changed.
- Jets: `SetJets(pt, eta, phi, mass, btagDeepFlavB, chEmEF, neEmEF, chHEF, neHEF,
  muEF, nConstituents, chMultiplicity, neMultiplicity)` &mdash; v15 has no stored
  jet ID / PU-ID, so the worker recomputes it in `H4LTools::PassJetIDv15(i, isPUPPI)`.
- `SetGenJets(pt, eta, phi, mass)` is new (PU-jet / genuine-jet studies).
- FatJet handling (and hence most of the 2l2q resolved/boosted path) is commented
  out in `analyze` on this branch.

### MET-phi correction

Done in Python inside `H4LCppModule.analyze`, not in a separate module. It
constructs `METPhiCorrector(campaign=Campaign.UL_20XX, is_data=..., is_puppi=True)`
from the tools' `met_phi_correction`, applies it to `PuppiMET.pt/phi` with
`npv=event.PV_npvs, run=event.run`, and the corrected values feed the worker
(`corr_pt`, `corr_phi`) and the `pT_MET` / `phi_MET` output branches. `Campaign`
only covers UL 2016/2017/2018 here &mdash; a `year` outside that set leaves
`corrector` undefined and raises.

### Config &mdash; `config/Input_<year>.yml`

Per-year: integrated `lumi`; `TriggerChannels` (list of keys, each a list of
`event.HLT_*` expressions evaluated with `eval` by
`modules/Helper.py:PassTrig`; each key also becomes a boolean output branch);
object cuts (`Electron`, `Muon`, `FsrPhoton`, `Jet` incl. deepJet b-tag WPs);
event cuts (`MZ1cut`, `MZZcut`, `Higgscut`, `MZcut`, `Zmass`); channel blocks
`HZZ2l2q` / `HZZ2l2nu`.

### Other modules

- `modules/Helper.py` &mdash; `PassTrig` + legacy python object-ID helpers.
- `modules/METFilters.py` &mdash; `passFilters(event, year)`.
- `modules/GenVarsProducer.py` &mdash; gen-level branches; flagged FIXME / "not
  working" in `post_proc.py`.
- `modules/JetSFMaker.py` &mdash; PU-ID scale factors (uses `data/`).
- `modules/keep_and_drop_list.py` &mdash; `keep_drop_rules_GEN` /
  `keep_drop_rules_Data_MC` lists.
- `modules/H4Lmodule.py` &mdash; older pure-python producer, **not** used by
  `post_proc.py`.

### Output

One hadd-ed skimmed NanoAOD file plus `cutFlow.json` and a `cutFlow` `TH1F`
written into the output file. Cut-flow counts accumulate on the C++ worker (one
attribute per `dynamicCuts_*` entry) and are dumped in `endFile`/`endJob`.

## Downstream: `higgs_combine/`

Standalone ROOT C++ macros (not part of the skim, not built by `scram`) for limit
extraction with the Higgs Combine tool. Run with ROOT, e.g.
`root -l -b -q higgs_combine/make_shapes.cpp`. `make_shapes.cpp` /
`make_shapes_final.cpp` build the shape-input ROOT file from skimmed ntuples
(cross sections, generated-event counts, categories hard-coded in the macros);
`make_datacard.cpp` reads that file's histogram yields and writes the Combine
datacard.

## Other scripts

`scripts/analysis/` &mdash; DY stitching, duplicate removal, JSON merge, sync text,
and lumi dumping (`dump_lumi.py` and `LumiDumper.py`, a `Module` that collects
run&rarr;lumiblock sets and writes JSON). `scripts/plotting/PlotCutFlowHistogram.py`.
`scripts/condor/`, `scripts/crab/`, `scripts/utils/`. Loose top-level utilities:
`scripts/fileList.py`, `scripts/check_das_sample.py`, `scripts/inspectNanoFile.py`,
`scripts/print_hist_content.py`.

## Conventions & gotchas

- `.gitignore` excludes `*.sh`, `*.txt`, `*.dat`, `*.root`, `*.jdl`, `*.patch`,
  `__init__.py`, `summary.dat`, `JHUGenMELA/`, `yaml-cpp/`. Generated condor
  files, most sample lists, and the two `external/` sub-repos are untracked;
  `git add -f` or edit `.gitignore` to commit one.
- `post_proc.py` has CRLF line endings; keep them if editing it.
- `SyncLepton2018GGH.txt` is opened for writing in `H4LCppModule.beginFile` &mdash;
  the CWD must be writable.
- `docs/README.md` still references some pre-reorg paths
  (`condor_setup_lxplus.py`, `scripts/GetLogSummary.py`) that now live under
  `scripts/`.
