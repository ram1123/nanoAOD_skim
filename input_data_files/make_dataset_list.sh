#!/bin/bash
year=${1:-2016}
fileout=sample_list_v7_${year}_eos_custom.dat

echo "# SM EWK" > ${fileout}
eos root://cmseos.fnal.gov find -d --maxdepth 3 /eos/uscms/store/group/lnujj/VVjj_aQGC/custom_nanoAOD | cut -d '=' -f 2 | cut -d '/' -f 1,9-11 | \
    grep -E "EWK_LO_SM.*${year}.*v7.*2010" >> ${fileout}
echo "# SM QCD" >> ${fileout}
eos root://cmseos.fnal.gov find -d --maxdepth 3 /eos/uscms/store/group/lnujj/VVjj_aQGC/custom_nanoAOD | cut -d '=' -f 2 | cut -d '/' -f 1,9-11 | \
    grep -E "QCD_LO_SM.*${year}.*v7.*2010" >> ${fileout}
echo "# aQGC" >> ${fileout}
eos root://cmseos.fnal.gov find -d --maxdepth 3 /eos/uscms/store/group/lnujj/VVjj_aQGC/custom_nanoAOD | cut -d '=' -f 2 | cut -d '/' -f 1,9-11 | \
    grep -E "EWK_LO_aQGC.*${year}.*v7.*2010" >> ${fileout}

