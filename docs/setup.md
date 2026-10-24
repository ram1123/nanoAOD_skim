# Setup

The package expects to sit at a fixed path inside a CMSSW release. Almost every
path is resolved relative to that directory or to `$CMSSW_BASE`.

## One-time installation

```bash
wget https://raw.githubusercontent.com/ram1123/nanoAOD_skim/HZZ_Analysis_2l2nu/setup.sh
SKIM_BRANCH=HZZ_Analysis_2l2nu /cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER -- bash setup.sh
# setup.sh alone checks out HZZ_Analysis_2l2nu_dev; SKIM_BRANCH picks the branch
```

[`setup.sh`](https://github.com/ram1123/nanoAOD_skim/blob/HZZ_Analysis_2l2nu/setup.sh)
does the following:

1. `scram project CMSSW_14_0_2`.
2. Clones the **`nanoAOD-tools` fork** (`ram1123/nanoAOD-tools`, branch
   `h4l_allCh_dev`). This fork provides `modules/common/met_phi_correction.py`
   (`METPhiCorrector`, `Campaign`), which upstream lacks.
3. Clones this package.
4. Builds `external/yaml-cpp` as a shared library (patched with
   `external/yamlcpp_pkg_py2to3.patch`).
5. Builds **JHUGenMELA** `v2.4.2` by hand: its own `setup.sh` calls `tcsh`, which
   the AF does not have. This step also downloads `libmcfm_710.so`,
   `libMG_SMEFTsim_v1.so` and the NNPDF grid from `spin.pha.jhu.edu`.
6. Copies the b-tag CSVs to `PhysicsTools/NanoAODTools/data/btagSF/`.
7. Stashes an `openssl-1.1` shim in `external/el8compat_lib/`.
8. Runs `scram b`.

`src/H4LTools.cc` is **not** compiled by `scram`. ROOT ACLiC JIT-compiles it on
first use (`.L H4LTools.cc+O`), so the first run after a C++ change is slower.

## Why the `cmssw-el9` container

`CMSSW_14_0_2` is an el8 build, but the pre-compiled JHU `libmcfm_710.so` needs
glibc &ge; 2.29. el8 hosts (the Purdue AF, Hammer, lxplus8) have glibc 2.28, so
on those hosts the build and every run happen inside `cmssw-el9`:

```bash
/cvmfs/cms.cern.ch/common/cmssw-el9 -B /depot/cms/users/$USER -B /work/users/$USER -- bash run.sh
```

| gotcha | fix |
| --- | --- |
| The wrapper silently drops bind paths that are not in the base image, so the script is "No such file" (exit 127) | pass every area explicitly with `-B` |
| `bash -c '<multi-line>'` mangles newlines | put commands in a script file |
| The container banner leaks into scram's environment probe, so `cmsenv` silently fails (`CMSSW_BASE` empty) | `export SCRAMRT_SET=1` before `scramv1 runtime -sh` |
| el8-built ROOT needs `libssl.so.1.1`, which AlmaLinux 9 lacks | prepend `external/el8compat_lib` to `LD_LIBRARY_PATH` before `import ROOT` |

On lxplus9 the container is not needed.

## Environment script for a run

```bash
# run.sh: execute with  cmssw-el9 -B /depot/cms/users/$USER -B /work/users/$USER -- bash run.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh; export SCRAM_ARCH=el8_amd64_gcc12
export SCRAMRT_SET=1
cd /work/users/$USER/<area>/CMSSW_14_0_2/src && eval `scramv1 runtime -sh`
cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
export X509_USER_PROXY=$HOME/x509_proxy          # voms-proxy-init -voms cms beforehand
export LD_LIBRARY_PATH=$PWD/external/el8compat_lib:$LD_LIBRARY_PATH
eval $(external/JHUGenMELA/MELA/setup.sh env)
python3 post_proc.py -i config/ExampleInputFileList.txt -n 300
```

`set_env.sh` is the lxplus recipe (MELA `LD_LIBRARY_PATH` plus `voms-proxy-init`).
It does not enter a container, so it will not work standalone on the AF.
