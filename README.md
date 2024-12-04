# NanoAOD Skim
nanoAOD skiming code for H->ZZ->4l/2l2Q/2l2nu studies.

# Setup

To setup the code, download the setup script and run it. As it will download four GitHub repository using ssh link, so it will ask for the password. The details of the `setup.sh` script is given in this [README.md](docs/README.md) file.

```bash
wget https://raw.githubusercontent.com/ram1123/nanoAOD_skim/refs/heads/HZZ_Analysis_Merge/setup.sh
source setup.sh
```

Once it's done, you can set the environment and run the code.

```bash
cmsenv
source set_env.sh
python3 post_proc.py
```

# batch job submission.
Condor-job submission (recommended)

In the file [condor_setup_lxplus.py](scripts/condor/condor_setup_lxplus.py), specify the correct input text file (present inside directory [input_data_Files](input_data_Files)) from which you need to take input NanoAOD DAS names. Also, updated the output EOS path. Then do the following:

   ```bash
   cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim
   # Use the arguments that you need.
   python3 scripts/condor/condor_setup_lxplus.py --input_file sample_list_v9_2018.dat
   # Set proxy before submitting the condor jobs.
   voms-proxy-init -voms cms --valid 200:00
   condor_submit <Files-created-from-above-command>.jdl
   ```

   To resubmit the failed jobs, use the following command:

   ```bash
   python3 scripts/condor/nanoAOD_condor_resubmit.py -d condor_logs/SkimNanoAOD_2022_ZXCR/240312_135155/ -s /eos/user/r/rasharma/nanoAOD_ntuples/SkimNanoAOD_2022_ZXCR/ -i submit_condor_jobs_lnujj_SkimNanoAOD_2022_ZXCR.txt -n 1
   ```

   This will give you new txt file. Then you can submit the condor job using new txt file.
