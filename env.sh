#!/bin/bash

# Update LD_LIBRARY_PATH for JHUGenMELA
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/JHUGenMELA/MELA/data/el9_amd64_gcc12

# Initialize a new proxy with the desired validity
voms-proxy-init --voms cms --valid 168:00

if [ $? -eq 0 ]; then
    echo "Proxy successfully created."

    # Check if the proxy is created in /tmp
    PROXY_PATH=$(voms-proxy-info --path)

    if [[ $PROXY_PATH == /tmp/* ]]; then
        echo "Proxy is located in /tmp, moving it to home directory..."
        echo "cp $PROXY_PATH ~/"
        cp $PROXY_PATH ~/
        echo "export X509_USER_PROXY=~/$(basename $PROXY_PATH)"
        export X509_USER_PROXY=~/$(basename $PROXY_PATH)
        echo "Proxy moved to home directory and X509_USER_PROXY set to $X509_USER_PROXY"
    else
        echo "Proxy is not in /tmp, no need to move it."
    fi
else
    echo "Failed to create the proxy."
    exit 1
fi
