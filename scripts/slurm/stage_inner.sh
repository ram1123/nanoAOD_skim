#!/bin/bash
# Runs INSIDE cmssw-el9. Mirror the CMSSW work area to /depot (Slurm cannot see
# /work) and pre-build the ROOT ACLiC lib so array tasks only load it.
set -uo pipefail

SRC=${1:?src CMSSW}
DST=${2:?dst CMSSW}
PKGREL=src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim

echo "### rsync $SRC/  ->  $DST/"
mkdir -p "$DST"
rsync -a --delete \
  --exclude='.git/' --exclude='*.pyc' --exclude='__pycache__/' \
  --exclude='tmp/' --exclude='logs/' --exclude='Pdfdata' \
  --exclude='src/H4LTools_cc*' \
  "$SRC/" "$DST/" || { echo "rsync failed"; exit 1; }

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12

# 1) relocate the release's internal paths BEFORE reading its runtime env
cd "$DST/src"
scram b ProjectRename >/dev/null 2>&1 && echo "### ProjectRename ok" || echo "### ProjectRename returned nonzero (continuing)"

# 2) now load the runtime env -- CMSSW_BASE must point at the /depot copy
eval "$(scramv1 runtime -sh 2>/dev/null)"
echo "### CMSSW_BASE=$CMSSW_BASE"
if [ "$CMSSW_BASE" != "$DST" ]; then
  echo "### ERROR: CMSSW_BASE ($CMSSW_BASE) != staged path ($DST)"
  echo "###        the copied release did not relocate; aborting."
  exit 1
fi

PKG="$DST/$PKGREL"
cd "$PKG"
export X509_USER_PROXY="${X509_USER_PROXY:-$HOME/x509_proxy}"
export LD_LIBRARY_PATH="$PKG/external/el8compat_lib:${LD_LIBRARY_PATH:-}"
eval "$(external/JHUGenMELA/MELA/setup.sh env)"

echo "### pre-build src/H4LTools_cc.so"
rm -f src/H4LTools_cc* src/H4LTools_cc.d
python3 - "$DST" <<'PY'
import ROOT, sys
dst = sys.argv[1]
from modules.H4LCppModule import HZZAnalysisCppProducer
HZZAnalysisCppProducer(year=2018, cfgFile="config/Input_2018.yml", isMC=True,
                       isFSR=True, cutFlowJSONFile="pb.json", channels="2l2v", DEBUG=False)
so = [l for l in ROOT.gSystem.GetLibraries().split() if "H4LTools_cc" in l]
print("### loaded from:", so)
assert so and so[0].startswith(dst), "ACLiC lib not under the staged /depot path: %s" % so
PY
rc=$?
rm -f pb.json SyncLepton2018GGH.txt
echo "### resulting src/:" ; ls -la src/H4LTools_cc.so 2>&1
[ $rc -eq 0 ] && [ -f src/H4LTools_cc.so ] && echo "### PREBUILD OK -> $DST" || { echo "### PREBUILD FAILED"; exit 1; }

echo "### tar the staged release for fast node-local extraction in jobs"
TGZ="$(dirname "$DST")/$(basename "$DST").tgz"
tar --exclude=.git -C "$(dirname "$DST")" -czf "$TGZ.tmp" "$(basename "$DST")" && mv -f "$TGZ.tmp" "$TGZ"
ls -la "$TGZ"
echo "### staged -> $DST"
echo "### tarball -> $TGZ  (jobs extract this to node-local /tmp)"
