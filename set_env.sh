#!/bin/bash
# ------------------------------------------------------------------------------
# Per-session environment for the HZZ NanoAOD skim.
#
#   SOURCE this file from anywhere, *inside* a cmssw-el9 shell:
#
#       /cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER --interactive   # (AF / lxplus8)
#       source <path>/nanoAOD_skim/set_env.sh
#       python3 post_proc.py -i config/ExampleInputFileList.txt
#
#   It does:  cmsenv  +  MELA env  +  openssl-1.1 shim  +  VOMS proxy.
#   On a real el9 host (lxplus9) the container line is not needed.
# ------------------------------------------------------------------------------

# --- locate the package + release from this script's own path ---------------
_SE_SRC="${BASH_SOURCE[0]:-$0}"
PKG="$(cd "$(dirname "$_SE_SRC")" && pwd)"
# .../CMSSW_x/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
RELEASE="$(cd "$PKG/../../../../../../.." && pwd)"

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH="${SCRAM_ARCH:-el8_amd64_gcc12}"

# --- cmsenv --------------------------------------------------------------
# `cmsenv` / bare `scramv1 runtime -sh` are broken inside the Purdue-AF
# cmssw-el9 image: its welcome banner leaks into scram's internal
# `env -i bash -c 'command -v python3'` and scram exits 127. Pre-setting
# SCRAMRT_SET to non-empty makes scram take the clean `which python3` path.
if [ -z "${CMSSW_BASE:-}" ] || [ "$CMSSW_BASE" != "$RELEASE" ]; then
    echo "set_env.sh: cmsenv ($RELEASE)"
    ( cd "$RELEASE/src" ) || { echo "set_env.sh: $RELEASE/src not found" >&2; return 1 2>/dev/null || exit 1; }
    pushd "$RELEASE/src" >/dev/null
    export SCRAMRT_SET="${SCRAMRT_SET:-1}"
    eval "$(scramv1 runtime -sh 2>/dev/null)"
    popd >/dev/null
fi
if [ -z "${CMSSW_BASE:-}" ]; then
    echo "set_env.sh: cmsenv FAILED" >&2
    return 1 2>/dev/null || exit 1
fi

# --- MELA -------------------------------------------------------------------
# sets MELA_ARCH / MELA_LIB_PATH / LD_LIBRARY_PATH / PYTHONPATH / ROOFITSYS
if [ -x "$PKG/external/JHUGenMELA/MELA/setup.sh" ]; then
    echo "set_env.sh: MELA environment"
    eval "$(cd "$PKG" && external/JHUGenMELA/MELA/setup.sh env)"
else
    echo "set_env.sh: WARNING - $PKG/external/JHUGenMELA/MELA missing; run setup.sh" >&2
fi

# --- openssl 1.1 shim -----------------------------------------------------
# el8-built ROOT needs libssl.so.1.1 / libcrypto.so.1.1; the cmssw-el9
# container (AlmaLinux 9) ships v3. setup.sh stashes a host copy here.
if [ -d "$PKG/external/el8compat_lib" ]; then
    export LD_LIBRARY_PATH="$PKG/external/el8compat_lib:${LD_LIBRARY_PATH:-}"
fi

# --- grid proxy --------------------------------------------------------
if voms-proxy-info -exists --valid 1:00 >/dev/null 2>&1; then
    echo "set_env.sh: VOMS proxy OK"
else
    echo "set_env.sh: voms-proxy-init"
    voms-proxy-init --rfc --voms cms --valid 168:00
fi
_PROXY="/tmp/x509up_u$(id -u)"
if [ -f "$_PROXY" ]; then
    cp -f "$_PROXY" "$HOME/x509_proxy" 2>/dev/null && export X509_USER_PROXY="$HOME/x509_proxy"
    echo "set_env.sh: X509_USER_PROXY=$X509_USER_PROXY"
else
    echo "set_env.sh: WARNING - no proxy at $_PROXY" >&2
fi

echo "set_env.sh: ready.  cd $PKG ;  python3 post_proc.py -i config/ExampleInputFileList.txt"
