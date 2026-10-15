#!/usr/bin/env python3
"""Merge skimmed NanoAOD files with haddnano.py.

    python3 scripts/analysis/mergeNanoAODRootFiles.py -i <inputDir> -o <outputDir> -f <output.root> [-r]
    python3 scripts/analysis/mergeNanoAODRootFiles.py -i skims/<submission> -o <outputDir> --sample-dirs
        (every sub-folder with ROOT files -> <outputDir>/<sub-folder>.root)

Needs ROOT (run after cmsenv, inside cmssw-el9 on the Purdue AF).

- Every input is checked (opens, not zombie/recovered, has an Events and a Runs tree).
  Unusable inputs are listed and the merge stops with exit code 1, unless --allow-skip
  (skipping a data file silently loses luminosity).
- Any haddnano.py failure stops the script with a non-zero exit code.
- The output is rebuilt unless it exists, is valid, and was made from exactly the same
  input files (recorded in <output>.inputs.txt); --force always rebuilds.
- More than --batch-size inputs are merged in batches inside a fresh temporary directory
  that is removed afterwards; the output is written to <output>.part and renamed only
  after a successful merge.
"""
import argparse
import glob
import logging
import os
import shutil
import subprocess
import sys
import tempfile

import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))


def find_haddnano(user_path=None):
    """haddnano.py from --haddnano, PATH (after cmsenv), $CMSSW_BASE, or this checkout."""
    candidates = [user_path, shutil.which("haddnano.py")]
    if os.environ.get("CMSSW_BASE"):
        candidates.append(os.path.join(os.environ["CMSSW_BASE"], "src/PhysicsTools/NanoAODTools/scripts/haddnano.py"))
    # this file: .../NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/scripts/analysis
    candidates.append(os.path.join(HERE, *[".."] * 6, "scripts", "haddnano.py"))
    for c in candidates:
        if c and os.path.isfile(c):
            return os.path.abspath(c)
    sys.exit("ERROR: haddnano.py not found (run cmsenv, or pass --haddnano <path>)")


def check_file(fname):
    """Return None if the file is usable for merging, otherwise the reason it is not."""
    if not os.path.exists(fname):
        return "does not exist"
    try:
        f = ROOT.TFile.Open(fname)  # recent PyROOT raises OSError instead of returning null
    except OSError as e:
        return "cannot be opened ({})".format(str(e).splitlines()[0])
    if not f:
        return "cannot be opened"
    try:
        if f.IsZombie():
            return "zombie"
        if f.TestBit(ROOT.TFile.kRecovered):
            return "recovered (file was not closed properly)"
        for tree in ("Events", "Runs"):
            obj = f.Get(tree)
            if not obj or not obj.InheritsFrom("TTree"):
                return "no {} tree".format(tree)
        # an Events tree with 0 entries is fine: the Runs tree (genEventSumw) must still be merged
        return None
    finally:
        f.Close()


def run_haddnano(haddnano, target, inputs):
    cmd = [haddnano, target] + list(inputs)
    logging.info("haddnano.py %s  <- %d files", target, len(inputs))
    logging.debug("Command: %s", " ".join(cmd))
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        raise RuntimeError("haddnano.py failed (exit code {}) while writing {}".format(rc, target))
    problem = check_file(target)
    if problem:
        raise RuntimeError("haddnano.py output {} is not usable: {}".format(target, problem))


def merge_files(haddnano, target, inputs, batch_size, tmp_base):
    """Merge into <target>.part, then rename to target. Temporary files live in a fresh
    directory that is always removed."""
    part = target + ".part"
    if os.path.exists(part):
        os.remove(part)
    if len(inputs) <= batch_size:
        run_haddnano(haddnano, part, inputs)
    else:
        tmpdir = tempfile.mkdtemp(prefix="mergeNano_", dir=tmp_base)
        logging.info("%d files: merging in batches of %d in %s", len(inputs), batch_size, tmpdir)
        try:
            batches = []
            for i in range(0, len(inputs), batch_size):
                batch_target = os.path.join(tmpdir, "batch{}.root".format(len(batches)))
                run_haddnano(haddnano, batch_target, inputs[i:i + batch_size])
                batches.append(batch_target)
            run_haddnano(haddnano, part, batches)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    os.replace(part, target)


class MergeError(Exception):
    pass


def merge_one(input_dir, target, args, haddnano, recursive):
    """Merge the ROOT files of input_dir into target. Returns "merged" or "up to date";
    raises MergeError on any problem."""
    record = target + ".inputs.txt"
    if recursive:
        found = glob.glob(os.path.join(input_dir, "**", "*.root"), recursive=True)
    else:
        found = glob.glob(os.path.join(input_dir, "*.root"))
    # never merge the output (or its partial file) into itself; sorted for a reproducible order
    inputs = sorted(os.path.abspath(f) for f in found if os.path.abspath(f) not in (target, target + ".part"))
    if not inputs:
        raise MergeError("no ROOT files found in {}".format(input_dir))
    logging.info("Found %d ROOT files in %s", len(inputs), input_dir)

    bad = []
    for f in inputs:
        problem = check_file(f)
        if problem:
            bad.append((f, problem))
    if bad:
        for f, problem in bad:
            logging.error("Unusable input: %s (%s)", f, problem)
        if not args.allow_skip:
            raise MergeError("{} of {} inputs unusable; fix or resubmit them, or rerun with --allow-skip".format(
                len(bad), len(inputs)))
        logging.warning("--allow-skip: merging %d of %d files", len(inputs) - len(bad), len(inputs))
        bad_files = set(f for f, _ in bad)
        inputs = [f for f in inputs if f not in bad_files]
        if not inputs:
            raise MergeError("no usable inputs left")

    if not args.force and os.path.exists(target) and os.path.exists(record):
        with open(record) as fh:
            previous = [line.strip() for line in fh if line.strip()]
        if previous == inputs and check_file(target) is None:
            logging.info("Up to date (same %d inputs): %s -- use --force to rebuild", len(inputs), target)
            return "up to date"
        logging.info("Inputs changed since the last merge: rebuilding %s", target)

    try:
        merge_files(haddnano, target, inputs, args.batch_size, args.tmpdir)
    except RuntimeError as e:
        raise MergeError(str(e))
    with open(record, "w") as fh:
        fh.write("\n".join(inputs) + "\n")
    logging.info("Merged %d files into %s (input list: %s)", len(inputs), target, record)
    return "merged"


def main():
    parser = argparse.ArgumentParser(
        description="Merge skimmed NanoAOD ROOT files with haddnano.py.",
        epilog="Single folder:  -i <dir> -o <outDir> -f <out.root>.   "
               "All samples:    -i skims/<submission> -o <outDir> --sample-dirs  "
               "(every sub-folder with ROOT files -> <outDir>/<sub-folder>.root)")
    parser.add_argument("-i", "--inputDir", required=True,
                        help="Folder with the ROOT files to merge (with --sample-dirs: folder of sample folders).")
    parser.add_argument("-o", "--outputDir", required=True, help="Directory for the merged file(s).")
    parser.add_argument("-f", "--outputFile", help="Name of the merged ROOT file (not used with --sample-dirs).")
    parser.add_argument("-r", "--recursive", action="store_true", help="Also search sub-directories (single-folder mode).")
    parser.add_argument("--sample-dirs", action="store_true",
                        help="Merge each sub-folder of --inputDir into <outputDir>/<sub-folder>.root.")
    parser.add_argument("--force", action="store_true", help="Rebuild outputs even if they are up to date.")
    parser.add_argument("--allow-skip", action="store_true",
                        help="Merge the usable files even if some inputs are unusable (default: stop).")
    parser.add_argument("--batch-size", type=int, default=100, help="Max files per haddnano.py call.")
    parser.add_argument("--tmpdir", default=os.environ.get("TMPDIR", "/tmp"),
                        help="Base directory for batch temp files (a fresh sub-directory is used).")
    parser.add_argument("--haddnano", default=None, help="Path to haddnano.py (default: auto-detect).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(levelname)s: %(message)s")

    if not os.path.isdir(args.inputDir):
        sys.exit("ERROR: input directory does not exist: {}".format(args.inputDir))
    if not args.sample_dirs and not args.outputFile:
        sys.exit("ERROR: -f/--outputFile is required (or use --sample-dirs)")
    os.makedirs(args.outputDir, exist_ok=True)
    haddnano = find_haddnano(args.haddnano)
    logging.info("haddnano.py: %s", haddnano)

    if not args.sample_dirs:
        target = os.path.abspath(os.path.join(args.outputDir, args.outputFile))
        try:
            merge_one(args.inputDir, target, args, haddnano, args.recursive)
        except MergeError as e:
            sys.exit("ERROR: {}".format(e))
        return 0

    # --sample-dirs: one merged file per sub-folder that contains ROOT files
    samples = sorted(d for d in os.listdir(args.inputDir)
                     if os.path.isdir(os.path.join(args.inputDir, d))
                     and glob.glob(os.path.join(args.inputDir, d, "*.root")))
    if not samples:
        sys.exit("ERROR: no sub-folders with ROOT files in {}".format(args.inputDir))
    logging.info("%d sample folders in %s", len(samples), args.inputDir)
    status = {}
    for n, sample in enumerate(samples, 1):
        logging.info("===== [%d/%d] %s", n, len(samples), sample)
        target = os.path.abspath(os.path.join(args.outputDir, sample + ".root"))
        try:
            status[sample] = merge_one(os.path.join(args.inputDir, sample), target, args, haddnano, False)
        except MergeError as e:
            logging.error("%s: %s", sample, e)
            status[sample] = "FAILED: {}".format(e)

    failed = [s for s, st in status.items() if st.startswith("FAILED")]
    logging.info("===== Summary: %d merged, %d up to date, %d failed",
                 sum(st == "merged" for st in status.values()),
                 sum(st == "up to date" for st in status.values()), len(failed))
    for s in failed:
        logging.error("  %s -> %s", s, status[s])
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
