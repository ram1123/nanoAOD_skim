#!/bin/bash

_set_env_fail() {
    local message="$1"
    echo "ERROR: $message" >&2
    return 1 2>/dev/null || exit 1
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MELA_LIB_DIR="$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/JHUGenMELA/MELA/data/el9_amd64_gcc12"
MELA_SETUP_SCRIPT="$SCRIPT_DIR/external/JHUGenMELA/MELA/setup.sh"

# Set LD_LIBRARY_PATH for MELA
echo "Updating LD_LIBRARY_PATH for MELA..."
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}$MELA_LIB_DIR"

if [ ! -f "$MELA_SETUP_SCRIPT" ]; then
    _set_env_fail "MELA setup script not found at $MELA_SETUP_SCRIPT"
fi

eval "$("$MELA_SETUP_SCRIPT" env)" || _set_env_fail "Failed to configure the MELA environment"

# Initialize VOMS proxy
echo "Initializing VOMS proxy..."
voms-proxy-init --voms cms --valid 168:00 || _set_env_fail "voms-proxy-init failed"

# Copy the proxy file to the home directory
PROXY_FILE=/tmp/x509up_u$(id -u)
if [ -f "$PROXY_FILE" ]; then
    echo "Copying proxy file to the home directory..."
    echo "Proxy file: $PROXY_FILE"
    echo "Home directory: $HOME"
    echo "Copying proxy file to $HOME/x509_proxy"
    cp "$PROXY_FILE" $HOME/x509_proxy
    export X509_USER_PROXY=$HOME/x509_proxy
    echo "Proxy file copied and X509_USER_PROXY is set."
else
    echo "Proxy file not found! Please ensure voms-proxy-init was successful."
    _set_env_fail "Proxy file not found after voms-proxy-init"
fi

echo "Environment setup complete. Proxy is active, and all settings are configured."
