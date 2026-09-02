#!/bin/bash
# Runs INSIDE cmssw-el9 on a Hammer node. One NanoAOD file -> one skim.
# Env in (from job_wrapper.sh): HZZ_CMSSW HZZ_LFN HZZ_SAMPLE HZZ_OUTDIR HZZ_TAG
#                              HZZ_PROXY [HZZ_ENTRIES] [HZZ_CHANNELS] [HZZ_SYST]
set -euo pipefail

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12

CMSSW=${HZZ_CMSSW:?}
PKG="$CMSSW/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"

cd "$CMSSW/src"
eval "$(scramv1 runtime -sh 2>/dev/null)"
cd "$PKG"

export X509_USER_PROXY=${HZZ_PROXY:?}
# el8-built ROOT needs openssl 1.1, absent on AlmaLinux 9
export LD_LIBRARY_PATH="$PKG/external/el8compat_lib:${LD_LIBRARY_PATH:-}"
eval "$(external/JHUGenMELA/MELA/setup.sh env)"

WORK="/tmp/hzzskim_${SLURM_JOB_ID:-x}_${SLURM_ARRAY_TASK_ID:-0}"
rm -rf "$WORK"; mkdir -p "$WORK/external"; cd "$WORK"
# post_proc.py + H4LTools.cc open these by path *relative to CWD*
ln -sfn "$PKG/config"                        config
ln -sfn "$PKG/data"                          data
ln -sfn "$PKG/external/CoupleConstantsForMELA" external/CoupleConstantsForMELA

OUT="${HZZ_SAMPLE}_${HZZ_TAG}_Skim.root"
CF="cutFlow_${HZZ_SAMPLE}_${HZZ_TAG}.json"
ENTRIES=${HZZ_ENTRIES:-0}
CHANNELS=${HZZ_CHANNELS:-2l2v}
SYST=""; [ -n "${HZZ_SYST:-}" ] && SYST="--WithSyst"

echo "=== $(date) task $HZZ_TAG  sample=$HZZ_SAMPLE ==="
echo "    in : $HZZ_LFN"
echo "    out: $HZZ_OUTDIR/$OUT"

set -x
python3 "$PKG/post_proc.py" -i "$HZZ_LFN" -outDir "$WORK" -o "$OUT" -c "$CF" \
        --entriesToRun "$ENTRIES" --channels "$CHANNELS" $SYST
set +x

test -s "$WORK/$OUT" || { echo "ERROR: $OUT not produced"; exit 2; }
mkdir -p "$HZZ_OUTDIR"
cp -v "$WORK/$OUT" "$HZZ_OUTDIR/$OUT"
cp -v "$WORK/$CF"  "$HZZ_OUTDIR/$CF" || true

cd /tmp; rm -rf "$WORK"
echo "=== $(date) DONE task $HZZ_TAG ==="
