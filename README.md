# NanoAOD Skim
nanoAOD skiming code for H->ZZ->2l2Q studies.

## Code setup

1. Step: 1: Get CMSSW release

   ```bash
   cmsrel CMSSW_14_0_2
   cd CMSSW_14_0_2/src
   cmsenv
   ```

2. Step: 2: Get  official nanoAODTools

   ```bash
   git clone git@github.com:cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
   cd PhysicsTools/NanoAODTools
   git checkout d163c18096fe2c5963ff5a9764bb420b46632178 # Updated to commit on 6 Dec 2023 in official nanoAOD-tools
   ```

3. Step: 3: Get our analysis repository

   ```bash
   cd $CMSSW_BASE/src
   git clone git@github.com:ram1123/nanoAOD_skim.git PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   cd PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   git checkout ZXCR
   cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools
   git apply python/postprocessing/analysis/nanoAOD_skim/external/nanoAODTools_py2to3.patch
   cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   git clone git@github.com:jbeder/yaml-cpp.git external/yaml-cpp
   cd external/yaml-cpp/
   git apply ../yamlcpp_pkg_py2to3.patch
   mkdir build
   cd build
   cmake3 .. -DBUILD_SHARED_LIBS=ON
   cmake3 --build .
   cd $CMSSW_BASE/src
   cp PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/btag/*.csv PhysicsTools/NanoAODTools/data/btagSF/.
   scram b
   voms-proxy-init --voms cms --valid 168:00
   ```

   (Optional: Fix git repo)

   ```bash
   find PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/.git/ -name "*.py*" -delete
   ```

4. Step: 4: Get the MELA package

   ```bash
   cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   git clone -b v2.4.2 https://github.com/JHUGen/JHUGenMELA
   cd JHUGenMELA
   git apply ../external/JHUGen_py2to3.patch
   cd ..
   sh JHUGenMELA/MELA/setup.sh -j 8
   cd JHUGenMELA/MELA/data/el9_amd64_gcc12/
   chmod +x *.so
   ```

4. Step: 4: interactive running

   ```bash
   cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/r/rasharma/h2l2nu/checkNewSetup_15July2024/CMSSW_14_0_2/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/JHUGenMELA/MELA/data/el9_amd64_gcc12
   # NOTE: The above export command is needed to run just before running the post_proc.py script. Otherwise, it will give error.
   python3 post_proc.py
   ```

5. batch job submission.
   1. Step: 5 (a): Condor-job submission (recommended)
      1. In the file [condor_setup_lxplus.py](condor_setup_lxplus.py), specify the correct input text file (present inside directory [input_data_Files](input_data_Files)) from which you need to take input NanoAOD DAS names. Also, updated the output EOS path. Then do the following:

         ```bash
         cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
         # Use the arguments that you need.
         python condor_setup_lxplus.py --input_file sample_list_v12_2022.dat
         # Set proxy before submitting the condor jobs.
         voms-proxy-init -voms cms --valid 200:00
         condor_submit <Files-created-from-above-command>.jdl
         ```

         To resubmit the failed jobs, use the following command:

         ```bash
         python Utils/nanoAOD_condor_resubmit.py -d condor_logs/SkimNanoAOD_2022_ZXCR/240312_135155/ -s /eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCR/ -i submit_condor_jobs_lnujj_SkimNanoAOD_2022_ZXCR.jdl -n 1
         ```

         This will give you new jdl file. Then you can submit the new jdl file.

   1. Step: 5(b): Crab-job submission (Not tested recently)
      ```bash
      cd crab/
      voms-proxy-init -voms cms --valid 200:00
      source /cvmfs/cms.cern.ch/crab3/crab.sh
      crab submit -c crab_cfg.py
      ```

## Few additioanl scripts

1. [condor_setup_lxplus.py](condor_setup_lxplus.py): This script can be used to setup the condor jobs. It takes the input text file (present inside directory [input_data_Files](input_data_Files)) from which you need to take input NanoAOD DAS names. Also, updated the output EOS path. Then do the following:

   ```bash
   python condor_setup_lxplus.py --input-file sample_list_v12_2022.dat

   # or
   python condor_setup_lxplus.py --submission_name SkimNanoAOD_2022_ZXCR --input_file sample_list_v12_2022.dat --condor_queue tomorrow
   ```

   This will create the condor job files and the condor log files.

1. [scripts/GetLogSummary.py](scripts/GetLogSummary.py): This script can be used to get the summary of the condor jobs. It takes the condor log files as input and gives the summary of the jobs. This summary contains the cut-flow table. It can be used as follows:

   ```bash
   python scripts/GetLogSummary.py <condor_log_file_base_path>
   ```

2. [scripts/check_das_sample.py](scripts/check_das_sample.py): This script can be used to check the status of the DAS samples. It takes the DAS name of the sample as input and gives the status of the sample. It can be used as follows:

   ```bash
   python scripts/check_das_sample.py <DAS_name_of_the_sample>
   ```

3. [scripts/nanoAOD_condor_resubmit.py](scripts/nanoAOD_condor_resubmit.py): This script can be used to resubmit the failed condor jobs. It takes the condor log files as input and resubmits the failed jobs. It can be used as follows:

   ```bash
   python nanoAOD_condor_resubmit.py -d <condor_log_file_base_path> -s <output_eos_path> -i <submit_jdl_file> -n <number_of_jobs_to_submit>

   # Example command:
   python nanoAOD_condor_resubmit.py -d condor_logs/SkimNanoAOD_2022_v12/240229_091018 -s /eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_v12/ -i submit_condor_jobs_lnujj_SkimNanoAOD_2022_v12.jdl -n 1
   ```

4. [scripts/mergeNanoAODRootFiles.py](scripts/mergeNanoAODRootFiles.py): This script can be used to merge the nanoAOD root files. It takes the input directory and the output directory as input and merges the nanoAOD root files. It can be used as follows:

   ```bash
   python scripts/mergeNanoAODRootFiles.py -i <input_directory> -o <output_directory> -f <output_file_name>

   # Example command:
   python scripts/mergeOutput.py -i /eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCR/EGamma/Run2022G/240312_135155/ -o /eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCR/EGamma -f Run2022G.root
   ```


## Few important points
