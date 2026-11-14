#!/bin/bash
# ==============================================================================
# One-time build of the HZZ NanoAOD skim and RestFrames.
#
# Run it INSIDE an el9 environment, from an empty working directory:
#
#     wget https://raw.githubusercontent.com/ram1123/nanoAOD_skim/HZZ_Analysis_2l2nu_dev/setup.sh
#     /cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/$USER -- bash setup.sh
#
# On a real el9 host (lxplus9) you can drop the `cmssw-el9 ... --` wrapper.
# On lxplus8 / the Purdue Analysis Facility (el8 host) the wrapper is required:
# CMSSW_14_0_2 is an el8 build but the JHU pre-compiled MCFM library needs
# glibc >= 2.29, so everything runs in the el9 container.
#
# Override the checked-out branch with:  SKIM_BRANCH=<branch> bash setup.sh
# ==============================================================================
set -e

CMSSW_RELEASE=CMSSW_14_0_2
export SCRAM_ARCH=el8_amd64_gcc12
SKIM_BRANCH="${SKIM_BRANCH:-HZZ_Analysis_2l2nu}"
TOOLS_BRANCH="${TOOLS_BRANCH:-h4l_allCh_dev}"
MELA_TAG="${MELA_TAG:-v2.4.2}"
NJOBS="${NJOBS:-8}"
RESTFRAMES_COMMIT=6fa7fbe356258d48654c5bc0a76f61d55cd9201c

# use https clones (ssh keys are usually not set up on lxplus / the AF)
GH=https://github.com

echo "=============================================================="
echo " release      : $CMSSW_RELEASE  ($SCRAM_ARCH)"
echo " skim branch  : $SKIM_BRANCH"
echo " tools branch : $TOOLS_BRANCH"
echo " MELA tag     : $MELA_TAG"
echo " OS           : $(cat /etc/redhat-release 2>/dev/null),  glibc $(ldd --version | head -1 | grep -oE '[0-9.]+$')"
echo "=============================================================="

source /cvmfs/cms.cern.ch/cmsset_default.sh

# The Purdue-AF cmssw-el9 image prints a welcome banner that leaks into scram's
# internal `env -i bash -c 'command -v python3'` and makes bare `cmsenv` /
# `scramv1 runtime -sh` exit 127. Pre-setting SCRAMRT_SET dodges that path.
export SCRAMRT_SET="${SCRAMRT_SET:-1}"

# ---- 1. CMSSW release ----------------------------------------------------
echo ">>> [1/8] scram project $CMSSW_RELEASE"
[ -d "$CMSSW_RELEASE" ] || scram project CMSSW "$CMSSW_RELEASE"
cd "$CMSSW_RELEASE/src"
eval "$(scramv1 runtime -sh 2>/dev/null)"
[ -n "${CMSSW_BASE:-}" ] || { echo "cmsenv failed"; exit 1; }
echo "    CMSSW_BASE=$CMSSW_BASE"

# ---- 2. nanoAOD-tools (fork) ------------------------------------------
echo ">>> [2/8] nanoAOD-tools ($TOOLS_BRANCH)"
[ -d PhysicsTools/NanoAODTools ] || \
    git clone -b "$TOOLS_BRANCH" "$GH/ram1123/nanoAOD-tools.git" PhysicsTools/NanoAODTools

# ---- 3. this analysis package ----------------------------------------
PKGREL=PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
echo ">>> [3/8] nanoAOD_skim ($SKIM_BRANCH)"
[ -d "$PKGREL" ] || \
    git clone -b "$SKIM_BRANCH" "$GH/ram1123/nanoAOD_skim.git" "$PKGREL"
PKG="$CMSSW_BASE/src/$PKGREL"

# ---- 4. yaml-cpp -----------------------------------------------------
echo ">>> [4/8] yaml-cpp"
cd "$PKG"
if [ ! -f external/yaml-cpp/build/libyaml-cpp.so ]; then
    rm -rf external/yaml-cpp
    git clone "$GH/jbeder/yaml-cpp.git" external/yaml-cpp
    ( cd external/yaml-cpp && git apply ../yamlcpp_pkg_py2to3.patch || true
      mkdir -p build && cd build
      cmake .. -DBUILD_SHARED_LIBS=ON -DYAML_CPP_BUILD_TESTS=OFF
      cmake --build . -j "$NJOBS" )
fi

# ---- 5. b-tag CSVs -------------------------------------------------
echo ">>> [5/8] b-tag CSV files"
mkdir -p "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF"
cp -f data/btag/*.csv "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF/" || true

# ---- 6. MELA -------------------------------------------------------
# JHUGenMELA/MELA/setup.sh cannot be used here: its dependency step calls
# `tcsh` (not installed on the AF). Do the equivalent by hand.
echo ">>> [6/8] JHUGenMELA ($MELA_TAG)"
cd "$PKG/external"
if [ ! -f JHUGenMELA/MELA/data/el9_amd64_gcc12/libJHUGenMELAMELA.so ]; then
    [ -d JHUGenMELA ] || git clone -b "$MELA_TAG" "$GH/JHUGen/JHUGenMELA"
    ( cd JHUGenMELA && git apply ../JHUGen_py2to3.patch || true )

    cd JHUGenMELA/MELA
    MELA_ARCH=el9_amd64_gcc12                 # MELA's own name for gcc12, any OS
    export MELA_ARCH
    export MELA_LIB_PATH="$PWD/data/$MELA_ARCH"
    export ROOFITSYS="${ROOFITSYS:-$ROOTSYS}"
    export LD_LIBRARY_PATH="$MELA_LIB_PATH:${LD_LIBRARY_PATH:-}"
    mkdir -p "$MELA_LIB_PATH"

    # (a) COLLIER  (wget + make, no tcsh)
    [ -f "$MELA_LIB_PATH/libcollier.so" ] || sh COLLIER/setup.sh -j "$NJOBS"

    # (b) NNPDF grid  (bash script)
    [ -f data/Pdfdata/NNPDF30_lo_as_0130.LHgrid ] || bash downloadNNPDF.sh

    # (c) pre-compiled MCFM + MadGraph SMEFTsim  (what tcsh retrieve.csh would fetch)
    MCFM_URL="$(cat "$MELA_LIB_PATH/download.url")"
    for lib in libmcfm_710.so libMG_SMEFTsim_v1.so; do
        [ -f "$MELA_LIB_PATH/$lib" ] || \
            wget --no-check-certificate -q -O "$MELA_LIB_PATH/$lib" "${MCFM_URL}${lib}"
    done

    # (d) build JHUGen (fortran) then the MELA C++ wrapper
    ( cd fortran && make -j "$NJOBS" )        # -> libjhugenmela.so
    make nopython -j "$NJOBS"                 # -> libJHUGenMELAMELA.so
    chmod +x "$MELA_LIB_PATH"/*.so
fi

# ---- 6b. openssl 1.1 shim for the el9 container --------------------
# el8-built ROOT links libssl.so.1.1 / libcrypto.so.1.1; AlmaLinux 9 ships v3.
# Stash the el8 host copies so set_env.sh can put them on LD_LIBRARY_PATH.
SHIM="$PKG/external/el8compat_lib"
if [ ! -e "$SHIM/libssl.so.1.1" ]; then
    mkdir -p "$SHIM"
    for base in libssl libcrypto; do
        src=$(ls /usr/lib64/${base}.so.1.1* 2>/dev/null | head -1)
        if [ -n "$src" ]; then
            cp -f "$src" "$SHIM/" && ln -sf "$(basename "$src")" "$SHIM/${base}.so.1.1"
        fi
    done
    [ -e "$SHIM/libssl.so.1.1" ] && echo "    openssl-1.1 shim -> $SHIM" \
        || echo "    (host has no openssl 1.1 -- fine on a real el9 host)"
fi

# ---- 7. scram build ----------------------------------------------
echo ">>> [7/8] scram b -j $NJOBS"
cd "$CMSSW_BASE/src"
scram b -j "$NJOBS"

# ---- 8. RestFrames (same ROOT/compiler as CMSSW; no separate LCG setup) ----
echo ">>> [8/8] RestFrames ($RESTFRAMES_COMMIT)"
cd "$PKG"
RF_SRC="$PKG/external/RestFrames-source"
export RESTFRAMESSYS="$PKG/external/RestFrames"
RF_CXX="$(command -v c++)"
RF_ROOT_CONFIG="$(command -v root-config)"
command -v rootcling
case "$RF_CXX" in /cvmfs/cms.cern.ch/*) ;; *)
    echo "ERROR: compiler is outside the CMS environment: $RF_CXX" >&2; exit 1;; esac
case "$RF_ROOT_CONFIG" in /cvmfs/cms.cern.ch/*) ;; *)
    echo "ERROR: ROOT is outside the CMS environment: $RF_ROOT_CONFIG" >&2; exit 1;; esac
export LD_LIBRARY_PATH="$SHIM:${LD_LIBRARY_PATH:-}"
if [ -e "$RESTFRAMESSYS" ]; then
    [ -f "$RESTFRAMESSYS/lib/libRestFrames.so" ] || {
        echo "ERROR: incomplete RestFrames installation: $RESTFRAMESSYS" >&2; exit 1;
    }
    echo "    Reusing existing RestFrames; testing loading, not certifying build provenance."
else
    [ ! -e "$RF_SRC" ] || {
        echo "ERROR: $RF_SRC already exists; inspect it before retrying a clean build." >&2
        exit 1
    }
    git clone -b dev/XZllZnn "$GH/ram1123/RestFrames.git" "$RF_SRC"
    git -C "$RF_SRC" checkout --detach "$RESTFRAMES_COMMIT"
    ( cd "$RF_SRC"
      CXX="$RF_CXX" bash ./configure --prefix="$RESTFRAMESSYS" --disable-debug
      make -j "$NJOBS"
      make install )
    [ -f "$RESTFRAMESSYS/lib/libRestFrames.so" ] || exit 1
    printf 'Commit: %s\nROOT: %s\nCompiler: %s\n' \
        "$RESTFRAMES_COMMIT" "$(root-config --version)" "$RF_CXX" \
        > "$RESTFRAMESSYS/build_environment.txt"
fi
export LD_LIBRARY_PATH="$RESTFRAMESSYS/lib:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$RESTFRAMESSYS/include:${ROOT_INCLUDE_PATH:-}"
python3 - <<'PY'
import os
import ROOT
lib = os.path.join(os.environ['RESTFRAMESSYS'], 'lib', 'libRestFrames.so')
if ROOT.gSystem.Load(lib) < 0:
    raise RuntimeError('Cannot load RestFrames: ' + lib)
ROOT.gInterpreter.AddIncludePath(os.path.join(os.environ['RESTFRAMESSYS'], 'include'))
if not ROOT.gInterpreter.Declare(r'''
#include "RestFrames/RestFrames.hh"
void check_rjr_installation() {
    RestFrames::ContraBoostInvJigsaw cb("install_cb", "install_cb");
}
'''):
    raise RuntimeError('RestFrames headers could not be declared')
ROOT.check_rjr_installation()
print('PASS: RestFrames library, headers and ContraBoostInvJigsaw')
PY

# Generate a relocatable environment entry point inside the installed package.
# Refuse to overwrite an existing different script.
RF_ENV_TMP="$(mktemp "$PKG/.set_env_rjr.XXXXXX")"
cat > "$RF_ENV_TMP" <<'RJR_ENV'
#!/bin/bash
# Source this from the installed package, inside the framework EL9 environment.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "Use: source set_env_rjr.sh" >&2
    exit 1
fi
_RJR_PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
export SCRAM_ARCH=el8_amd64_gcc12
source "$_RJR_PKG/set_env.sh" || return 1
export RESTFRAMESSYS="$_RJR_PKG/external/RestFrames"
if [ ! -f "$RESTFRAMESSYS/lib/libRestFrames.so" ]; then
    echo "Missing RestFrames library under $RESTFRAMESSYS" >&2
    return 1
fi
export PATH="$RESTFRAMESSYS/bin:$PATH"
export LD_LIBRARY_PATH="$RESTFRAMESSYS/lib:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$RESTFRAMESSYS/include:${ROOT_INCLUDE_PATH:-}"
echo "RestFrames environment ready: $RESTFRAMESSYS"
unset _RJR_PKG
RJR_ENV
if [ -e "$PKG/set_env_rjr.sh" ]; then
    if ! cmp -s "$RF_ENV_TMP" "$PKG/set_env_rjr.sh"; then
        echo "ERROR: existing set_env_rjr.sh differs. Review proposed file: $RF_ENV_TMP" >&2
        exit 1
    fi
    rm -f "$RF_ENV_TMP"
else
    mv "$RF_ENV_TMP" "$PKG/set_env_rjr.sh"
fi

cat <<EOF

==============================================================
 Build complete.

 Each session:
   /cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/\$USER -- bash   # AF / lxplus8
   source /cvmfs/cms.cern.ch/cmsset_default.sh
   export SCRAM_ARCH=el8_amd64_gcc12
   cd $CMSSW_BASE/src && cmsenv
   cd $PKGREL
   source set_env_rjr.sh
   python3 post_proc.py -i config/ExampleInputFileList.txt

 See README.md and scripts/slurm/README.md (batch on Purdue).
 RestFrames is installed under external/RestFrames.
 RJR reconstruction helper/module edits must be applied separately if absent.
==============================================================
EOF

