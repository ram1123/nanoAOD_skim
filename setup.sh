#!/bin/bash

# Setup CMSSW release
echo "Setting up CMSSW release..."
cmsrel CMSSW_14_0_2
cd CMSSW_14_0_2/src
cmsenv

# Clone official nanoAODTools repository
echo "Cloning official nanoAODTools repository..."
git clone -b h4l_allCh_dev git@github.com:ram1123/nanoAOD-tools.git PhysicsTools/NanoAODTools

# Clone and setup analysis repository
echo "Setting up analysis repository..."
git clone -b HZZ_Analysis_2l2nu git@github.com:ram1123/nanoAOD_skim.git PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim

# Build YAML-cpp dependency
echo "Setting up YAML-cpp..."
cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
git clone git@github.com:jbeder/yaml-cpp.git external/yaml-cpp
cd external/yaml-cpp/
git apply ../yamlcpp_pkg_py2to3.patch
mkdir build
cd build
cmake3 .. -DBUILD_SHARED_LIBS=ON
cmake3 --build .

# Copy BTag CSV files
echo "Copying BTag CSV files..."
cd $CMSSW_BASE/src
cp PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/btag/*.csv PhysicsTools/NanoAODTools/data/btagSF/.

# Clone and setup MELA package
echo "Setting up MELA package..."
cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/external
git clone -b v2.4.2 https://github.com/JHUGen/JHUGenMELA
cd JHUGenMELA
git apply ../JHUGen_py2to3.patch
cd ..
sh JHUGenMELA/MELA/setup.sh -j 8
cd JHUGenMELA/MELA/data/el9_amd64_gcc12/
chmod +x *.so

# Compile the project
cd $CMSSW_BASE/src
echo "Compiling CMSSW..."
scram b -j 8

# Back to the working directory
cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/r/rasharma/h2l2nu/checkNewSetup_15July2024/CMSSW_14_0_2/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/external/JHUGenMELA/MELA/data/el9_amd64_gcc12

# Note for interactive running
echo "Setup for interactive running completed!"
echo "Before running post_proc.py, execute:"
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/r/rasharma/h2l2nu/checkNewSetup_15July2024/CMSSW_14_0_2/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/external/JHUGenMELA/MELA/data/el9_amd64_gcc12'
echo ""
echo "OR set the path and proxy using the script: set_env.sh"
echo "source set_env.sh"
