# Slurm submission (Purdue AF / Hammer)

Runs the HZZ NanoAOD skim as a Slurm **job array** on the `hammer-nodes`
partition, one array task per input ROOT file.

## Why it is not like the condor scripts

Constraints discovered on Hammer (`scripts/condor/` targets lxplus HTCondor and
does **not** apply here):

| Constraint | Consequence |
|---|---|
| Jobs see `/cvmfs`, `/depot`, `/tmp` only — **not `/work` or `/eos`** | the CMSSW area must live on `/depot`; input via XRootD; output to `/depot` |
| Hammer nodes are el8 / glibc 2.28; the JHU `libmcfm_710.so` needs glibc ≥ 2.29 | every job runs inside the **`cmssw-el9`** container |
| el8-built ROOT needs `libssl.so.1.1`, absent on AlmaLinux 9 | `external/el8compat_lib/` (host openssl-1.1 copy) is prepended to `LD_LIBRARY_PATH` |
| `cmssw-el9` drops `APPTAINER_BINDPATH` entries whose mount point is not in the base image | `/depot` is bound with an explicit `-B` flag |
| MELA/MCFM first-time init is CPU-heavy on shared Hammer cores (~10–15 min) | this is a fixed per-job cost — give `--time` headroom and skim whole files per job, not tiny `--entriesToRun` |

## One-time setup (re-run after any code / MELA change)

```bash
bash scripts/slurm/stage_to_depot.sh
# -> rsyncs $CMSSW_BASE to /depot/cms/users/<you>/HZZ2l2nu_skim/CMSSW_14_0_2
#    runs `scram b ProjectRename`, pre-builds src/H4LTools_cc.so so array tasks
#    never recompile in the shared tree.
```

## Submit

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh        # for dasgoclient
voms-proxy-init -voms cms --valid 192:00           # proxy lifetime bounds the run

# dry run first: expand the DAS list, write job.sh + tasks.tsv, do NOT sbatch
python3 scripts/slurm/slurm_setup.py \
    --input_file sample_list_v15_2018.dat \
    --submission_name HZZ2l2nu_2018 \
    --max_files 2 --dry_run

# real submission
python3 scripts/slurm/slurm_setup.py \
    --input_file sample_list_v15_2018.dat \
    --submission_name HZZ2l2nu_2018
```

Key options (`--help` for all): `--output_base` (default
`/depot/cms/users/shar1172/HZZ2l2nu_skim`), `--cmssw_on_depot`, `--depot_bind`,
`--redirector` (default Purdue XCache), `--channels` (default `2l2v`),
`--entries` (0 = all), `--with_syst`, `--time`, `--mem`, `--max_parallel`
(array `%N` throttle), `--max_files` (debug cap).

## Layout produced

```
<output_base>/
  CMSSW_14_0_2/                        # staged release (from stage_to_depot.sh)
  submits/<name>_<timestamp>/
      job.sh          tasks.tsv        # the sbatch script + <idx>\t<lfn>\t<sample>\t<outdir>
      x509_proxy                       # copy Slurm can read
      logs/task_<n>.out                # per-task stdout+stderr
      summary.dat  git_*.patch         # provenance (best effort)
  skims/<name>/<sample>/<sample>_<idx>_Skim.root   (+ cutFlow_<sample>_<idx>.json)
```

`hadd` the per-file `*_Skim.root` per sample afterwards
(`scripts/analysis/mergeNanoAODRootFiles.py`).

## Files

| file | role |
|---|---|
| `stage_to_depot.sh` / `stage_inner.sh` | one-time /work → /depot mirror + H4LTools prebuild (runs in `cmssw-el9`) |
| `slurm_setup.py` | expand DAS list, generate `job.sh` + `tasks.tsv`, `sbatch` |
| `job_wrapper.sh` | array-task entry (host): pick tasks.tsv row, enter `cmssw-el9 -B <depot>` |
| `job_inner.sh` | inside container: cmsenv + MELA env + libssl shim, run `post_proc.py` in `/tmp`, copy skim to `/depot` |

## Resubmitting failures

`tasks.tsv` row index == array task id. Find failed tasks from the logs, then:

```bash
sbatch --array=<comma,list,of,ids> <output_base>/submits/<name>_<ts>/job.sh
```

## Known gaps

- Year/`isMC` detection in `post_proc.py` is substring-matching on the file path;
  `/store/data/RunXXXX/...` v15 paths may not carry a recognised token — check a
  data submission's first task log before scaling out.
- No automatic `hadd` / bookkeeping of processed vs total events yet.
