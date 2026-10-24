# Running the skim

## Interactive smoke test

There is no test suite. The canonical check is one run over the bundled one-file
list:

```bash
python3 post_proc.py -i config/ExampleInputFileList.txt -n 1000
```

`config/ExampleInputFileList.txt` currently holds a 2018A SingleMuon v15 data file
(1.67M events, about 57 min for the full file). Use a UL18 v15
`GluGluHToZZTo2L2Nu_M300` MC file to check MC weights and ee coverage.

## `post_proc.py` options

| option | default | meaning |
| --- | --- | --- |
| `-i/--inputFile` | `ExampleInputFileList.txt` in the **CWD** | a `.root` file, or a `.txt` with one LFN per line (`root://cms-xrd-global.cern.ch/` is prepended) |
| `-n/--entriesToRun` | `0` = all | event cap |
| `--channels` | `2l2v` | `2l2v`, `1g` (&gamma;+jets CR), `4l`, `2l2q`, `all` |
| `-o/--outputFile` | `skimmed_nano.root` | final hadd-ed file. Pass a **bare name**: it is joined with `-outDir` |
| `-outDir/--outputDir` | `.` | output directory |
| `-c/--cutFlowFile` | `cutFlow.json` | cut-flow dump (bare name) |
| `--noPreselection` | off | disables the ROOT-level lepton/photon preselection (slower, same output) |
| `--WithSyst` | off | adds the JME correctors. It **only stores** JES/JER-shifted branches; see [Corrections](corrections.md#not-applied-yet) |
| `--DEBUG` | off | per-event printout |

### Year and data/MC detection

The year is detected by **substring matching on the first input file's path**,
and the tokens are specific to the NanoAOD version. If nothing matches, the job
fails because `cfgFile` stays `None`.

| token in path | era | config | golden JSON |
| --- | --- | --- | --- |
| `UL16NanoAODAPVv15`, `HIPM_UL2016_NanoAODv15` | 2016preVFP | `Input_2016preVFP.yml` | `Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt` |
| `UL16NanoAODv15`, `UL2016_NanoAODv15` | 2016postVFP | `Input_2016postVFP.yml` | same |
| `UL17NanoAODv15`, `UL2017` | 2017 | `Input_2017.yml` | `Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt` |
| `UL18NanoAODv15`, `UL2018_NanoAODv15` | 2018 | `Input_2018.yml` | `Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt` |

`isMC = "/data/" not in path`.

## Batch: Slurm on Purdue (Hammer)

`scripts/slurm/` runs one array task per input file, inside `cmssw-el9`, from a
tarball unpacked on node-local `/tmp`. Slurm jobs cannot see `/work`, so the
release is mirrored to `/depot` first.

```bash
bash scripts/slurm/stage_to_depot.sh             # once, and after every code / MELA change
source /cvmfs/cms.cern.ch/cmsset_default.sh       # dasgoclient
voms-proxy-init -voms cms --valid 192:00          # the proxy lifetime bounds the run

python3 scripts/slurm/slurm_setup.py --input_file sample_list_AN2016_325_v15_2017.dat \
        --submission_name HZZ2l2nu_AN325_2017_MCdata --max_files 2 --dry_run   # check first
python3 scripts/slurm/slurm_setup.py --input_file sample_list_AN2016_325_v15_2017.dat \
        --submission_name HZZ2l2nu_AN325_2017_MCdata
# gamma+jets CR:
python3 scripts/slurm/slurm_setup.py --input_file sample_list_gammaCR_v15_2017.dat \
        --submission_name HZZ2l2nu_gammaCR_2017 --channels 1g
```

A bare `--input_file` name is looked up in `input_data_files/`. Outputs land in
`<output_base>/skims/<name>/<sample>/<sample>_<idx>_Skim.root` (plus
`cutFlow_<sample>_<idx>.json`). The default `<output_base>` is
`/depot/cms/users/shar1172/HZZ2l2nu_skim`. Each job spends a fixed ~10&ndash;15
min on MELA/MCFM initialisation, so skim whole files rather than small `-n`.
To resubmit, run `sbatch --array=<ids> <submit>/job.sh`.

Details: `scripts/slurm/README.md`.

## Batch: HTCondor on lxplus

```bash
python3 scripts/condor/condor_setup_lxplus.py \
        --input_file input_data_files/sample_list_v15_2018.dat \
        --submission_name HZZ2l2nu_<date> --condor_queue tomorrow \
        --condor_file_name submit_condor_HZZ2l2nu_<date>
voms-proxy-init -voms cms --valid 200:00
condor_submit <generated>.jdl
```

Here `--input_file` is taken relative to the CWD as-is. Resubmit with
`scripts/condor/nanoAOD_condor_resubmit.py`.

## Cut flow

Each job writes `cutFlow.json` and a `cutFlow` TH1F into the output. The general
rows are:

| row | meaning |
| --- | --- |
| `Total` | all input entries (within `-n`) |
| `PassPresel` | events that reached `analyze`: preselection, plus the golden JSON for data |
| `PassTrig` | events passing any trigger list |
| `PassMETFilters` | events passing the MET noise filters |
| `PassZZ2l2nuSelection` / `PassZZ2l2nu_emuCR_Selection` / `PassGammaCR_Selection` | events selected per region |

The per-region sections follow, each giving the % of the previous step: the
ee/&mu;&mu; SR, the e&mu; CR, the &gamma; CR, and the muon object cuts.
`Total` is counted before the golden JSON and the preselection. The
`PassTrig` and `PassMETFilters` rows count preselected events.
