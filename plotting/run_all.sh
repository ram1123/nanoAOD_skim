#!/bin/bash
# Reproduce the 2017 v15 control plots with the AN-2016/325 selection
# (low-MET CR, m_ll-sideband CR, emu CR, alpha-method sidebands).
# Runs on the Purdue AF session (el8 host) with ROOT from the shared pixi env; no CMSSW needed.
#   bash plotting/run_all.sh                     (from nanoAOD_skim; SKIM / BASE can be overridden from the env)
# Detached:  setsid nohup bash plotting/run_all.sh > plotting/logs/run_all.log 2>&1 < /dev/null &
set -e
cd "$(dirname "$0")"                          # nanoAOD_skim/plotting
PIXI=/work/pixi/global/.pixi/envs/default
PY=$PIXI/bin/python
ROOTEXE=$PIXI/bin/root
# pixi ROOT finds the C++ std headers via its own compiler (x86_64-conda-linux-gnu-c++) on PATH; without it
# cling/ACLiC fail with "cannot extract standard library include paths". A CMSSW env (cmsenv / set_env.sh)
# on LD_LIBRARY_PATH makes it abort with the cling "--optimize-regalloc" error. Clean both.
export PATH=$PIXI/bin:$PATH
unset LD_LIBRARY_PATH ROOTSYS ROOT_INCLUDE_PATH PYTHONPATH CMSSW_BASE
SKIM=${SKIM:-/depot/cms/users/shar1172/HZZ2l2nu_skim/skims/HZZ2l2nu_AN325_2017_MCdata}
BASE=${BASE:-/depot/cms/users/shar1172/HZZ2l2nu_skim/controlplots_2017_AN325}
mkdir -p logs

# 1. slim + merge per sample, MC sum(genEventSumw), data de-duplication -> $BASE/inputs
# $PY ../scripts/analysis/prepare_inputs.py --skimdir $SKIM --outdir $BASE/inputs --threads 8 > logs/prepare_inputs.log 2>&1

# 2. sample list (xsec, 1/sumw) + cut files + variable header
$PY make_config.py --inputs $BASE/inputs

# 3. compile once, then the regions in parallel
$ROOTEXE -l -b -q -e '.L myControlPlots.C+' > /dev/null
for r in LowMET mllSideband emu alphaSB; do
    $ROOTEXE -l -b -q "RunMacro_2017_v15.C(\"$r\", \"$BASE/plots\")" > logs/plots_$r.log 2>&1 &
done
wait

# 4. yield table (+ alpha from data)
$PY yield_table.py $BASE/plots > yields_2017.md

# 5. metadata next to the outputs: what produced them
bash write_metadata.sh "$SKIM" "$BASE/inputs" "$BASE/plots"
echo DONE
