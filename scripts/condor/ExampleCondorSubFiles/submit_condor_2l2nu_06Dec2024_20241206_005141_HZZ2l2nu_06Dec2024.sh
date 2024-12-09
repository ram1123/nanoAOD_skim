#!/bin/bash
echo "Starting job on " `date`
echo "Running on: `uname -a`"
echo "System software: `cat /etc/redhat-release`"
source /cvmfs/cms.cern.ch/cmsset_default.sh
echo "====> List input arguments : "
echo "1. nanoAOD ROOT file: ${1}"
echo "2. EOS path to store output root file: ${2}"
echo "3. EOS path from where we copy CMSSW: ${3}"
echo "4. Output root file name: ${4}"
echo "========================================="
echo "copy cmssw tar file from store area"
xrdcp -f root://eosuser.cern.ch/${3}/CMSSW_14_0_2.tgz .
tar -xf CMSSW_14_0_2.tgz
rm CMSSW_14_0_2.tgz
cd CMSSW_14_0_2/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/
rm *.root
scramv1 b ProjectRename
eval `scram runtime -sh`
echo "========================================="
echo "cat post_proc.py"
echo "..."
cat post_proc.py
echo "..."
echo "========================================="
output_file=${4}_hadd.root
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/external/JHUGenMELA/MELA/data/el9_amd64_gcc12
python3 post_proc.py --entriesToRun 0 --inputFile ${1} --outputFile ${output_file}  --DownloadFileToLocalThenRun True  --channel "2l2v"
echo "====> List root files : "
ls -ltrh *.root
ls -ltrh *.json
echo "====> copying *.root file to stores area..."
if ls ${output_file} 1> /dev/null 2>&1; then
    echo "File ${output_file} exists. Copy this."
    echo "xrdcp -f ${output_file} root://eosuser.cern.ch/${2}/${4}_Skim.root"
    xrdcp -f ${output_file} root://eosuser.cern.ch/${2}/${4}_Skim.root
    echo "xrdcp -f ${4}.json root://eosuser.cern.ch/${2}/cutFlow_${4}.json"
    xrdcp -f ${4}.json root://eosuser.cern.ch/${2}/cutFlow_${4}.json
else
    echo "Something wrong: file ${output_file} does not exists, please check the post_proc.py script."
fi
rm *.root
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_14_0_2
