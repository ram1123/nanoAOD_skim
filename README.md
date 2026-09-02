# NanoAOD Skim

NanoAOD skimming code for H&rarr;ZZ&rarr;4l / 2l2q / 2l2nu studies.
Current development branch `HZZ_Analysis_2l2nu_dev` targets the **2l2nu** channel
and **NanoAOD v15** input.

## Platform note

This code is built for **el9** (AlmaLinux 9 / lxplus9). `CMSSW_14_0_2` is an el8
build, but the JHU pre-compiled MCFM library inside MELA needs `glibc >= 2.29`, so
on an **el8 host** (lxplus8, the Purdue Analysis Facility, Hammer) everything &mdash;
build *and* run &mdash; must happen inside the `cmssw-el9` container:

```bash
/cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER -- bash <command>
#   or, for an interactive shell:
/cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER --interactive
```

`-B <dir>` binds a writable path into the container &mdash; add every area you
need (`/work/users/$USER`, and `/depot/cms/users/$USER` for batch). On lxplus9 the
container is not needed; skip the `cmssw-el9 ... --` prefix everywhere below.

## Setup (one time)

```bash
wget https://raw.githubusercontent.com/ram1123/nanoAOD_skim/HZZ_Analysis_2l2nu_dev/setup.sh
/cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER -- bash setup.sh
```

[setup.sh](setup.sh) does: `scram project CMSSW_14_0_2` &rarr; clone the
`nanoAOD-tools` fork (branch `h4l_allCh_dev`) and this package &rarr; build
`yaml-cpp` &rarr; copy the b-tag CSVs &rarr; build **MELA** by hand (its own
`setup.sh` calls `tcsh`, which is absent on the AF) incl. downloading the
pre-compiled MCFM / MadGraph libs and the NNPDF grid from `spin.pha.jhu.edu`
&rarr; stash an `openssl-1.1` shim for the container &rarr; `scram b`.
Override the branch with `SKIM_BRANCH=<branch> bash setup.sh`.
Full step-by-step / lxplus notes: [docs/README.md](docs/README.md).

## Run interactively

```bash
/cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER --interactive     # el8 host only

REL=<path-to>/CMSSW_14_0_2/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
source $REL/set_env.sh              # cmsenv + MELA env + openssl-1.1 shim + VOMS proxy
cd $REL
python3 post_proc.py -i config/ExampleInputFileList.txt
```

[set_env.sh](set_env.sh) is self-contained &mdash; source it from anywhere inside a
`cmssw-el9` shell and it runs `cmsenv` for you (bare `cmsenv` is broken in the
Purdue-AF container: its welcome banner corrupts `scram`'s python detection; the
script works around it by pre-setting `SCRAMRT_SET`). On lxplus9, drop the
`cmssw-el9` line.

### `post_proc.py`

`config/ExampleInputFileList.txt` is a single-file list (a UL18 **v15**
`GluGluHToZZTo2L2Nu_M300` signal file) &mdash; the year (2018) and MC/data flag are
auto-detected from the path. The code on this branch is **NanoAOD v15 only**
(`Electron_mvaIso_WP90`, the v15 `SetJets` / `SetPuppiMET` signatures); a v9 file
fails with `RuntimeError: Unknown branch Electron_mvaIso_WP90`.

Common options (`python3 post_proc.py --help` for all):

| flag | meaning |
| --- | --- |
| `-i, --inputFile` | a `.root` file, or a `.txt` with one LFN per line (empty &rarr; reads `ExampleInputFileList.txt`) |
| `-n, --entriesToRun` | events to process; `0` = all (default) |
| `--channels {all,4l,2l2q,2l2v}` | default `2l2v` |
| `--DEBUG` | per-event selection printout |
| `-o / -outDir / -c` | output file / dir / cut-flow JSON name |

Output: a hadd-ed skim `skimmed_nano.root` + `cutFlow.json` + a `cutFlow` TH1F.

## Batch submission

### Slurm &mdash; Purdue AF / Hammer (recommended here)

`scripts/slurm/` &mdash; see **[scripts/slurm/README.md](scripts/slurm/README.md)**.
One-time `bash scripts/slurm/stage_to_depot.sh` (mirrors the release to `/depot`,
which Slurm can see and `/work` it cannot), then:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
voms-proxy-init --rfc --voms cms --valid 192:00
python3 scripts/slurm/slurm_setup.py --input_file sample_list_v15_2018.dat \
        --submission_name HZZ2l2nu_2018 --max_files 2 --dry_run   # inspect, then drop --dry_run
```

**Known issue:** the `Mela()` constructor currently hangs (>1 h, 100 % CPU) on
Hammer nodes &mdash; it is ~10 s interactively. The scaffolding is complete and jobs
start correctly, but do not scale out until this is resolved (a lazy/optional MELA
init for the 2l2nu channel, which never uses it, would sidestep it).

### HTCondor &mdash; lxplus

```bash
cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
python3 scripts/condor/condor_setup_lxplus.py --input_file input_data_files/sample_list_v15_2018.dat
voms-proxy-init --rfc --voms cms --valid 200:00
condor_submit <generated>.jdl
```

Resubmit failures with `scripts/condor/nanoAOD_condor_resubmit.py`.

## Sample lists

`input_data_files/sample_list_v15_<year>.dat` (MC) and
`..._<year>_data.dat` (data) for `2016preVFP`, `2016postVFP`, `2017`, `2018`,
regenerated from the copperhead dataset YAML by
`python3 input_data_files/make_run2_lists_from_yaml.py`. One DAS dataset per
line; `#` comments; `skip_sample` entries kept but commented.
