#!/bin/bash
# Slurm array-task entry point (runs on the Hammer node, host OS). Picks this
# task's line out of tasks.tsv, then hands off to job_inner.sh inside cmssw-el9.
#
# Env in (exported by the generated job.sh):
#   HZZ_CMSSW HZZ_TASKS HZZ_PROXY HZZ_BIND [HZZ_ENTRIES] [HZZ_CHANNELS] [HZZ_SYST]
set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?must be run as a job array}"
idx=$SLURM_ARRAY_TASK_ID

line=$(sed -n "$((idx + 1))p" "${HZZ_TASKS:?}")
[ -n "$line" ] || { echo "no task at row $idx of $HZZ_TASKS"; exit 1; }
IFS=$'\t' read -r T_IDX T_LFN T_SAMPLE T_OUTDIR <<< "$line"

export HZZ_LFN="$T_LFN" HZZ_SAMPLE="$T_SAMPLE" HZZ_OUTDIR="$T_OUTDIR" HZZ_TAG="$T_IDX"

INNER="$HZZ_CMSSW/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/scripts/slurm/job_inner.sh"
# -B is required: the cmssw-el9 wrapper drops APPTAINER_BINDPATH entries whose
# mount point is absent from the base image (/depot is), so bind it explicitly.
cd "${HZZ_CMSSW:?}"
exec /cvmfs/cms.cern.ch/common/cmssw-el9 -B "${HZZ_BIND:?}" -- bash "$INNER"
