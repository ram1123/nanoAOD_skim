# NanoAOD Skim

NanoAOD skimming code for H&rarr;ZZ&rarr;4l / 2l2q / 2l2nu studies, built on the
`nanoAOD-tools` `PostProcessor`. The 2l2nu development branches
(`HZZ_Analysis_2l2nu_dev`, `HZZ_Analysis_2l2nu_dev_gammaCR`) target the **2l2nu**
channel and **NanoAOD v15** input.

**Scope:** Run 2 UL **2017 and 2018** (`config/Input_2017.yml`, `Input_2018.yml`).
2016 is detected but `config/Input_2016.yml` does not exist yet.

## Workflow at a glance

| step | where it runs | tool |
| --- | --- | --- |
| 1. skim NanoAOD | `cmssw-el9` container (el8 hosts) | `post_proc.py`; at scale [Slurm](#slurm--purdue-af--hammer-recommended-here) or [HTCondor](#htcondor--lxplus) |
| 2. slim inputs per sample, Σ genEventSumw, data de-duplication | AF host, pixi ROOT | [`prepare_inputs.py`](#prepare-analysis-inputs-prepare_inputspy) |
| 3. control plots + yields | AF host, pixi ROOT | [`plotting/`](#control-plots-and-limits) |
| 4. shapes, datacard, limits | ROOT + Combine | [`higgs_combine/`](#control-plots-and-limits) |

Steps 2&ndash;4 must **not** run in a shell where `cmsenv` / `set_env.sh` was sourced
(see [Prepare analysis inputs](#prepare-analysis-inputs-prepare_inputspy)).

## Platform note

This code is built for **el9** (AlmaLinux 9 / lxplus9). `CMSSW_14_0_2` is an el8
build, but the JHU pre-compiled MCFM library inside MELA needs `glibc >= 2.29`, so
on an **el8 host** (lxplus8, the Purdue Analysis Facility, Hammer) everything in
step 1 &mdash; build *and* run &mdash; must happen inside the `cmssw-el9` container:

```bash
/cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER -- bash <script>
#   or, for an interactive shell:
/cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER --interactive
```

`-B <dir>` binds a path into the container; add every area you need
(`/work/users/$USER`, and `/depot/cms/users/$USER` for batch). Without it the
wrapper silently drops the mount and your script is "No such file" (exit 127).
Put multi-line commands in a script file: `bash -c '<multi-line>'` mangles newlines.
On lxplus9 the container is not needed; skip the `cmssw-el9 ... --` prefix below.

## Setup (one time)

```bash
wget https://raw.githubusercontent.com/ram1123/nanoAOD_skim/HZZ_Analysis_2l2nu_dev/setup.sh
/cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER -- bash setup.sh
```

[setup.sh](setup.sh) does: `scram project CMSSW_14_0_2` &rarr; clone the
`nanoAOD-tools` fork (`ram1123/nanoAOD-tools`, branch `h4l_allCh_dev`, which
provides the MET-&phi; corrector) and this package &rarr; build `yaml-cpp` &rarr;
copy the b-tag CSVs &rarr; build **MELA** by hand (its own `setup.sh` calls `tcsh`,
which is absent on the AF) incl. downloading the pre-compiled MCFM / MadGraph libs
and the NNPDF grid from `spin.pha.jhu.edu` &rarr; stash an `openssl-1.1` shim for
the container &rarr; `scram b`. Override the branch with
`SKIM_BRANCH=<branch> bash setup.sh`. Full step-by-step / lxplus notes:
[docs/README.md](docs/README.md).

`src/H4LTools.cc` is not built by `scram`: ROOT ACLiC compiles it on first use.

## Run interactively

```bash
/cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER --interactive     # el8 host only

REL=<path-to>/CMSSW_14_0_2/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
source $REL/set_env.sh              # cmsenv + MELA env + openssl-1.1 shim + VOMS proxy
cd $REL
python3 post_proc.py -i config/ExampleInputFileList.txt -n 1000
```

[set_env.sh](set_env.sh) is self-contained: source it from anywhere inside a
`cmssw-el9` shell and it runs `cmsenv` for you (bare `cmsenv` is broken in the
Purdue-AF container: its welcome banner corrupts `scram`'s python detection; the
script works around it by pre-setting `SCRAMRT_SET`). It refuses to run on an el8
shell. On lxplus9, drop the `cmssw-el9` line.

There is no test suite; the smoke test is one run over
`config/ExampleInputFileList.txt`. It holds a **2018A SingleMuon v15 data** file
(1.67M events, ~1 h in full, so use `-n`). For MC weights / ee coverage use a UL18
v15 signal file, e.g. `GluGluHToZZTo2L2Nu_M300` (in the git history of that file).
Other `config/ExampleInputFileList_*.txt` cover 4l / data / Run 3 / WW.

The code is **NanoAOD v15 only**: a v9 file fails with
`RuntimeError: Unknown branch Electron_mvaIso_WP90`.

### `post_proc.py`

| flag | meaning |
| --- | --- |
| `-i, --inputFile` | a `.root` file, or a `.txt` with one LFN per line (prefixed with `root://cms-xrd-global.cern.ch/`; empty &rarr; `ExampleInputFileList.txt`) |
| `-n, --entriesToRun` | events to process; `0` = all (default) |
| `--channels {all,4l,2l2q,2l2v,1g}` | default `2l2v`; `1g` = single-photon (&gamma;+jets) control-region skim |
| `--noPreselection` | drop the ROOT pre-cut (&ge; 2 leptons above the YAML p<sub>T</sub> cut &minus; 5 GeV). Output events are identical; only the `PassTrig` / `PassMETFilters` cut-flow rows change |
| `--WithSyst` | JME correctors + PU-ID SF. **Broken for MC** (`NameError` on `muonScaleRes` / `gammaSF`); leave it off |
| `--DEBUG` | per-event selection printout |
| `-o / -outDir / -c` | output file / dir / cut-flow JSON. Pass `-o` and `-c` as **bare file names**: `-o` is joined with `-outDir` for the final hadd |

Output: a hadd-ed skim (`skimmed_nano.root` by default) + `cutFlow.json` + a
`cutFlow` TH1F inside the file.

**Year and data/MC are detected from the first input path**, so LFNs must keep the
official dataset naming. MC vs data: `/data/` in the path. Year tokens:
`UL2017` / `UL17NanoAODv15` &rarr; 2017; `UL2018_NanoAODv15` / `UL18NanoAODv15`
&rarr; 2018; `20UL16NanoAOD(APV)v9` &rarr; 2016 (config missing); `Run2022` /
`Summer22` &rarr; 2022 (stale). No match &rarr; the job fails.

Cut values and triggers live in `config/Input_<year>.yml`, not in the C++.
Only unprescaled HLT paths are allowed there.

## Batch submission

### Slurm &mdash; Purdue AF / Hammer (recommended here)

`scripts/slurm/` &mdash; full details in **[scripts/slurm/README.md](scripts/slurm/README.md)**.
Slurm cannot see `/work`, so `stage_to_depot.sh` mirrors the release to `/depot`
(re-run it after every code change, but **not** while a run is in progress: tasks
extract the tarball when they start). Run on the AF host, not in the container:

```bash
bash scripts/slurm/stage_to_depot.sh                 # after code changes
source /cvmfs/cms.cern.ch/cmsset_default.sh          # dasgoclient
voms-proxy-init --rfc --voms cms --valid 192:00
python3 scripts/slurm/slurm_setup.py --input_file input_data_files/sample_list_v15_2017.dat \
        --submission_name HZZ2l2nu_2017 --redirector root://xcache.cms.rcac.purdue.edu/ \
        --max_files 2 --dry_run                      # inspect tasks.tsv, then drop --max_files/--dry_run
```

One array task per input file (4 h, 4 GB, 50 in parallel by default). Logs and
`tasks.tsv` are in `/depot/.../HZZ2l2nu_skim/submits/<name>_<timestamp>/`; outputs
in `/depot/.../HZZ2l2nu_skim/skims/<name>/<sample>/`. The proxy is copied into the
submit directory at submission, so renewing your own proxy does not renew it.

**Bad nodes:** some Hammer nodes hang every task until the wall time, with a log
holding only the two `couldn't chdir` lines (seen on `hammer-g002` and
`hammer-f010`). Group failures by node before blaming the code or raising `--time`.
Resubmit only the failed task ids, away from those nodes:

```bash
sbatch --exclude=hammer-f010,hammer-g002 --time=02:00:00 --array=<ids> <submit_dir>/job.sh
```

When the queue is empty, check that every `tasks.tsv` row has its
`<sample>_<id>_Skim.root` before moving on.

### HTCondor &mdash; lxplus

```bash
cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
python3 scripts/condor/condor_setup_lxplus.py --input_file input_data_files/sample_list_v15_2018.dat
voms-proxy-init --rfc --voms cms --valid 200:00
condor_submit <generated>.jdl
```

`--input_file` is taken as-is relative to the CWD. Resubmit failures with
`scripts/condor/nanoAOD_condor_resubmit.py`. `scripts/crab/` is an older,
untested path.

## Prepare analysis inputs (`prepare_inputs.py`)

Run once the skim jobs have finished, before plotting or making shapes.
[`scripts/analysis/prepare_inputs.py`](scripts/analysis/prepare_inputs.py) reads
`<skimdir>/<sample>/*_Skim.root` and writes to `<outdir>/`:

- one slim `<sample>.root` per MC sample: only events passing the 2l2nu SR or e&mu; CR
  selection, only the branches used downstream;
- `sumw.json`: &Sigma; genEventSumw per MC sample, summed over the `Runs` tree of
  **all** skim files (use this for normalisation, never an event count);
- one data file (`--data_file`, default `data_2017.root`): all primary datasets
  chained in priority order (DoubleMuon, DoubleEG, MuonEG, SingleMuon,
  SingleElectron, SinglePhoton) and de-duplicated on (run, lumi, event), the same
  logic as `removeDuplicates.C`.

It runs on the AF host with ROOT from the global pixi environment, **not** in
`cmssw-el9`. Use a fresh terminal: after `cmsenv` / `set_env.sh`, CMSSW's ROOT
libraries are on `LD_LIBRARY_PATH` and the pixi ROOT aborts with
`cling (LLVM option parsing): ... --optimize-regalloc ... may only occur zero or one times!`.
The `env -i ...` prefix below makes the command work from any shell. Keep the pixi `bin/` on
`PATH`: ROOT finds its C++ headers through the conda compiler there, and without it fails with
`cannot extract standard library include paths`. `plotting/run_all.sh` sets this up itself.

```bash
mkdir -p plotting/logs
env -i HOME=$HOME USER=$USER PATH=/work/pixi/global/.pixi/envs/default/bin:/usr/bin:/bin \
  /work/pixi/global/.pixi/envs/default/bin/python scripts/analysis/prepare_inputs.py \
    --skimdir /depot/cms/users/$USER/HZZ2l2nu_skim/skims/<submission_name> \
    --outdir  /depot/cms/users/$USER/HZZ2l2nu_skim/controlplots_2017_AN325/inputs \
    --threads 8 > plotting/logs/prepare_inputs.log 2>&1 &
tail -f plotting/logs/prepare_inputs.log
```

- Existing outputs are skipped, so re-running resumes. Pass `--force` (or clear
  `<outdir>`) after a re-skim.
- `--only <sample> ... | data` processes a subset.
- `--mode gamma` builds the single-photon CR inputs from a `--channels 1g` skim
  (SinglePhoton data only); run it separately with its own `--skimdir` / `--outdir`.
- The data file is built single-threaded (the de-duplication is not thread safe)
  and is the slowest part.

## Control plots and limits

Same environment as above (AF host, pixi ROOT, fresh terminal).

**Control plots:** [`plotting/`](plotting/README.md) (2017, AN-2016/325 selection:
low-MET, m<sub>ll</sub>-sideband, e&mu;, &alpha; sidebands, &gamma; CR).

```bash
cd /work/users/$USER/<area>/CMSSW_14_0_2/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim   # no cmsenv
export SKIM=/depot/cms/users/$USER/HZZ2l2nu_skim/skims/<submission_name>       # input: the Slurm skim (<SKIM>/<sample>/*_Skim.root)
export BASE=/depot/cms/users/$USER/HZZ2l2nu_skim/controlplots_<year>_<tag>     # output: $BASE/inputs (prepare_inputs.py --outdir), $BASE/plots
mkdir -p plotting/logs
setsid nohup bash plotting/run_all.sh > plotting/logs/run_all.log 2>&1 < /dev/null &
```

Set both every time. If unset, `run_all.sh` falls back to the last 2017 run
(`skims/HZZ2l2nu_AN325_2017_MCdata`, `controlplots_2017_AN325`), and a new skim
would silently be plotted from the old inputs. `BASE/inputs` must be the
`--outdir` you gave `prepare_inputs.py`.

`run_all.sh` runs, in order:

1. `prepare_inputs.py` (previous section). Samples already done are skipped, so this
   is quick if you ran it by hand first.
2. `plotting/make_config.py`: cross sections, cuts and ranges. It generates the
   sample list, `cuts/` and `controlplotvars_2017.h` (git-ignored). Edit these
   things here, not in the generated files.
3. `RunMacro_2017_v15.C` &rarr; `myControlPlots.C` (copied from PlottingCodes), then
   `yield_table.py` and `write_metadata.sh` (provenance next to the outputs on `/depot`).

**Limits:** [`higgs_combine/`](higgs_combine/) holds standalone ROOT macros
(`make_shapes*.cpp` &rarr; `make_datacard.cpp`). They currently read the full skims:
normalise MC by &Sigma; genEventSumw from the `Runs` tree, never by an event count
(genWeight can be ~100), and de-duplicate data across primary datasets first with
`root -l -b -q 'scripts/analysis/removeDuplicates.C+("inputs.txt","out.root")'`.

## Sample lists

`input_data_files/`, one DAS dataset per line, `#` = comment:

- `sample_list_v15_<year>.dat` (MC) and `..._<year>_data.dat` (data) for
  `2016preVFP`, `2016postVFP`, `2017`, `2018` (made from the copperhead dataset
  YAML; `skip_sample` entries are kept but commented).
- `sample_list_AN2016_325_v15_<year>[_data].dat`: the AN-2016/325 sample set
  (from `make_AN2016_325_lists.py`); `..._v9_...` are the NanoAOD v9 equivalents.
- `sample_list_gammaCR_v15_2017.dat`: inputs for the `--channels 1g` skim.

`.gitignore` excludes `*.dat`, `*.txt` and `*.sh` except the tracked ones; use
`git add -f` for a new list or script.

## Repository layout

| path | contents |
| --- | --- |
| `post_proc.py` | entry point: builds the module list and runs `PostProcessor` |
| `modules/H4LCppModule.py` | Python `Module`: reads the YAML, feeds collections to the C++ worker, writes branches and the cut flow |
| `src/H4LTools.cc`, `include/H4LTools.h` | object ID, FSR, ZZ / 2l2nu / e&mu; / &gamma; selections, MELA |
| `config/` | per-year YAML (cuts, triggers, lumi) and example input lists |
| `scripts/slurm/`, `scripts/condor/`, `scripts/crab/` | batch submission |
| `scripts/analysis/` | `prepare_inputs.py`, duplicate removal, DY stitching, lumi / JSON tools |
| `plotting/` | control-plot pipeline |
| `higgs_combine/` | shape inputs and datacards for Combine |
| `docs/` | full setup notes, AN summary, TODO |
