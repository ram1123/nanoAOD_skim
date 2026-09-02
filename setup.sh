#!/bin/bash
# ==============================================================================
# One-time build of the HZZ NanoAOD skim.
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
SKIM_BRANCH="${SKIM_BRANCH:-HZZ_Analysis_2l2nu_dev}"
TOOLS_BRANCH="${TOOLS_BRANCH:-h4l_allCh_dev}"
MELA_TAG="${MELA_TAG:-v2.4.2}"
NJOBS="${NJOBS:-8}"

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
echo ">>> [1/7] scram project $CMSSW_RELEASE"
[ -d "$CMSSW_RELEASE" ] || scram project CMSSW "$CMSSW_RELEASE"
cd "$CMSSW_RELEASE/src"
eval "$(scramv1 runtime -sh 2>/dev/null)"
[ -n "${CMSSW_BASE:-}" ] || { echo "cmsenv failed"; exit 1; }
echo "    CMSSW_BASE=$CMSSW_BASE"

# ---- 2. nanoAOD-tools (fork) ------------------------------------------
echo ">>> [2/7] nanoAOD-tools ($TOOLS_BRANCH)"
[ -d PhysicsTools/NanoAODTools ] || \
    git clone -b "$TOOLS_BRANCH" "$GH/ram1123/nanoAOD-tools.git" PhysicsTools/NanoAODTools

# ---- 3. this analysis package ----------------------------------------
PKGREL=PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
echo ">>> [3/7] nanoAOD_skim ($SKIM_BRANCH)"
[ -d "$PKGREL" ] || \
    git clone -b "$SKIM_BRANCH" "$GH/ram1123/nanoAOD_skim.git" "$PKGREL"
PKG="$CMSSW_BASE/src/$PKGREL"

# ---- 4. yaml-cpp -----------------------------------------------------
echo ">>> [4/7] yaml-cpp"
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
echo ">>> [5/7] b-tag CSV files"
mkdir -p "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF"
cp -f data/btag/*.csv "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF/" || true

# ---- 6. MELA -------------------------------------------------------
# JHUGenMELA/MELA/setup.sh cannot be used here: its dependency step calls
# `tcsh` (not installed on the AF). Do the equivalent by hand.
echo ">>> [6/7] JHUGenMELA ($MELA_TAG)"
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
echo ">>> [7/7] scram b -j $NJOBS"
cd "$CMSSW_BASE/src"
scram b -j "$NJOBS"

cat <<EOF

==============================================================
 Build complete.

 Each session:
   /cvmfs/cms.cern.ch/common/cmssw-el9 -B /work/users/\$USER -- bash   # AF / lxplus8
   source /cvmfs/cms.cern.ch/cmsset_default.sh
   export SCRAM_ARCH=el8_amd64_gcc12
   cd $CMSSW_BASE/src && cmsenv
   cd $PKGREL
   source set_env.sh
   python3 post_proc.py -i config/ExampleInputFileList.txt

 See README.md and scripts/slurm/README.md (batch on Purdue).
==============================================================
EOF
