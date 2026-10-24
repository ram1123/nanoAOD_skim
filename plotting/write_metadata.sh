#!/bin/bash
# Write <dir>/_metadata/ into the plotting inputs and plots folders: the configuration that made them
# (make_config.py, cuts, sample list, variable header, macros, yields), the git state of nanoAOD_skim,
# sum(genEventSumw), logs, and a copy of the skim's own _metadata (slurm_setup.py).
#   bash plotting/write_metadata.sh <skim_dir> <inputs_dir> <plots_dir>
set -e
SKIM=$1; shift
HERE=$(cd "$(dirname "$0")" && pwd)               # nanoAOD_skim/plotting
SK=$(cd "$HERE/.." && pwd)                        # nanoAOD_skim
for OUT in "$@"; do
    [ -d "$OUT" ] || continue
    M=$OUT/_metadata
    rm -rf "$M"; mkdir -p "$M"
    cp -r "$HERE/make_config.py" "$HERE/RunMacro_2017_v15.C" "$HERE/yield_table.py" "$HERE/run_all.sh" \
          "$HERE/myControlPlots.C" "$HERE/controlplotvars_2017.h" "$HERE/cuts" \
          "$SK/scripts/analysis/prepare_inputs.py" "$M/"
    cp "$HERE"/HZZ2l2nu_Samples13TeV_2017_v15*.txt "$M/" 2>/dev/null || true
    cp "$HERE/yields_2017.md" "$M/" 2>/dev/null || true
    cp -r "$HERE/logs" "$M/logs" 2>/dev/null || true
    [ -f "$(dirname "$OUT")/inputs/sumw.json" ] && cp "$(dirname "$OUT")/inputs/sumw.json" "$M/"
    [ -d "$SKIM/_metadata" ] && cp -r "$SKIM/_metadata" "$M/skim_metadata"
    {
        echo "created: $(date -Iseconds)  user: $USER  host: $(hostname)"
        echo "output : $OUT"
        echo "skims  : $SKIM"
        echo; echo "== nanoAOD_skim ($SK)"
        git -C "$SK" rev-parse --abbrev-ref HEAD; git -C "$SK" log -1 --format='%H %ad %s'
        git -C "$SK" status --short
    } > "$M/provenance.txt"
    git -C "$SK" diff HEAD > "$M/nanoAOD_skim_diff_HEAD.patch" || true
    echo "metadata -> $M"
done
