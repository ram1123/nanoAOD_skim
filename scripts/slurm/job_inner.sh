#!/bin/bash
# Runs INSIDE cmssw-el9 on a Hammer node. One NanoAOD file -> one skim.
#
# The staged /depot CMSSW tree is extracted to node-local /tmp first: MELA/MCFM
# init does thousands of small file ops and is ~50x slower over the /depot NFS
# than on local disk (10+ min vs seconds).
#
# Env in (from job_wrapper.sh): HZZ_CMSSW_TGZ HZZ_LFN HZZ_SAMPLE HZZ_OUTDIR
#                               HZZ_TAG HZZ_PROXY [HZZ_ENTRIES] [HZZ_CHANNELS] [HZZ_SYST]
set -uo pipefail

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
export PYTHONUNBUFFERED=1

JOBTMP="/tmp/hzzskim_${SLURM_JOB_ID:-x}_${SLURM_ARRAY_TASK_ID:-0}"
rm -rf "$JOBTMP"; mkdir -p "$JOBTMP"
trap 'cd /tmp; rm -rf "$JOBTMP"' EXIT

echo "=== $(date) task $HZZ_TAG  sample=$HZZ_SAMPLE  node=$(hostname) ==="
echo "    extracting $(basename "${HZZ_CMSSW_TGZ:?}") to $JOBTMP ..."
t0=$SECONDS
tar -xzf "$HZZ_CMSSW_TGZ" -C "$JOBTMP" || { echo "ERROR: tar extract failed"; exit 1; }
echo "    extracted in $((SECONDS-t0))s"

CMSSW="$JOBTMP/$(basename "${HZZ_CMSSW_TGZ%.tgz}")"
PKG="$CMSSW/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
[ -d "$PKG" ] || { echo "ERROR: $PKG missing after extract"; exit 1; }

cd "$CMSSW/src"
eval "$(scramv1 runtime -sh 2>/dev/null)"
[ -n "${CMSSW_BASE:-}" ] || { echo "ERROR: cmsenv failed"; exit 1; }
scram b ProjectRename >/dev/null 2>&1 || true
eval "$(scramv1 runtime -sh 2>/dev/null)"

cd "$PKG"
export X509_USER_PROXY=${HZZ_PROXY:?}
export LD_LIBRARY_PATH="$PKG/external/el8compat_lib:${LD_LIBRARY_PATH:-}"
eval "$(external/JHUGenMELA/MELA/setup.sh env)"

OUT="${HZZ_SAMPLE}_${HZZ_TAG}_Skim.root"
CF="cutFlow_${HZZ_SAMPLE}_${HZZ_TAG}.json"
ENTRIES=${HZZ_ENTRIES:-0}
CHANNELS=${HZZ_CHANNELS:-2l2v}
SYST=""; [ -n "${HZZ_SYST:-}" ] && SYST="--WithSyst"

echo "    in : $HZZ_LFN"
echo "    out: $HZZ_OUTDIR/$OUT"
set -x
python3 "$PKG/post_proc.py" -i "$HZZ_LFN" -outDir "$PKG" -o "$OUT" -c "$CF" \
        --entriesToRun "$ENTRIES" --channels "$CHANNELS" $SYST
rc=$?
set +x
[ $rc -eq 0 ] && [ -s "$PKG/$OUT" ] || { echo "ERROR: post_proc failed (rc=$rc) or no output"; exit 2; }

mkdir -p "$HZZ_OUTDIR"
cp -v "$PKG/$OUT" "$HZZ_OUTDIR/$OUT"
cp -v "$PKG/$CF"  "$HZZ_OUTDIR/$CF" || true
echo "=== $(date) DONE task $HZZ_TAG ==="
