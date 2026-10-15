#!/bin/bash
# One-time (re-run after any code / MELA change): mirror the CMSSW area to /depot,
# because Slurm jobs on Hammer cannot see /work.
#
#   bash scripts/slurm/stage_to_depot.sh
#   bash scripts/slurm/stage_to_depot.sh <SRC_CMSSW> <DST_CMSSW>
#
set -euo pipefail

SRC=${1:-${CMSSW_BASE:-/work/users/shar1172/XZZ2l2nu/CMSSW_14_0_2}}
DST=${2:-/depot/cms/users/shar1172/HZZ2l2nu_skim/CMSSW_14_0_2}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

[ -d "$SRC/src" ] || { echo "ERROR: $SRC is not a CMSSW release"; exit 1; }
SRC=$(cd "$SRC" && pwd)

# writable host subtrees that must appear inside the container (parents of SRC/DST)
SRC_BIND=/$(echo "$SRC"  | cut -d/ -f2-4)     # e.g. /work/users/shar1172
DST_BIND=/$(echo "$DST"  | cut -d/ -f2-5)     # e.g. /depot/cms/users/shar1172
echo "SRC (work) : $SRC   [bind $SRC_BIND]"
echo "DST (depot): $DST   [bind $DST_BIND]"

exec /cvmfs/cms.cern.ch/common/cmssw-el9 -B "$SRC_BIND" -B "$DST_BIND" -- \
     bash "$HERE/stage_inner.sh" "$SRC" "$DST"
