#!/usr/bin/env python3
"""
Generate (and optionally submit) a Slurm array job that runs the HZZ NanoAOD skim
on the Purdue AF / Hammer cluster.

Model
-----
* Slurm jobs on Hammer cannot see /work or /eos -- only /cvmfs, /depot, /tmp.
  So the CMSSW area must first be mirrored to /depot:

      bash scripts/slurm/stage_to_depot.sh

* Each dataset in the input list is expanded with `dasgoclient` to its files.
  One array task == one input ROOT file, read over XRootD (Purdue XCache by
  default), skimmed with post_proc.py inside a cmssw-el9 container, output copied
  to  <output_base>/skims/<submission>/<sample>/<sample>_<idx>_Skim.root .

Usage
-----
    source /cvmfs/cms.cern.ch/cmsset_default.sh   # (or cmsenv) -- for dasgoclient
    voms-proxy-init -voms cms --valid 192:00
    python3 scripts/slurm/slurm_setup.py --input_file sample_list_v15_2018.dat \
            --submission_name HZZ2l2nu_test --max_files 2 --dry_run
    # inspect, then:
    python3 scripts/slurm/slurm_setup.py --input_file sample_list_v15_2018.dat \
            --submission_name HZZ2l2nu_2018
"""
import argparse
import datetime
import os
import shutil
import subprocess
import sys

PKG_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# fixed location of this package inside any CMSSW release
PKG_REL = "src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
DASGOCLIENT = shutil.which("dasgoclient") or "/cvmfs/cms.cern.ch/common/dasgoclient"


def find_proxy():
    for p in (os.environ.get("X509_USER_PROXY"),
              "/tmp/x509up_u%d" % os.getuid(),
              os.path.expanduser("~/x509_proxy")):
        if p and os.path.isfile(p):
            return p
    return None


def das_files(dataset):
    out = subprocess.run([DASGOCLIENT, "--query", "file dataset=%s" % dataset],
                         capture_output=True, text=True)
    if out.returncode != 0:
        sys.stderr.write("dasgoclient failed for %s:\n%s\n" % (dataset, out.stderr))
        return []
    return [l.strip() for l in out.stdout.splitlines() if l.strip().endswith(".root")]


def sample_tag(dataset):
    # /Primary/Campaign-Cond/TIER  ->  Primary(+ext)
    parts = dataset.strip("/").split("/")
    name = parts[0] if parts else "sample"
    if "ext" in dataset:
        name += "_ext" + dataset.split("ext")[-1].split("/")[0]
    return name


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--input_file", required=True,
                    help="DAS-name list; bare name is looked up in input_data_files/")
    ap.add_argument("--submission_name", default="SkimNanoAOD")
    ap.add_argument("--output_base", default="/depot/cms/users/shar1172/HZZ2l2nu_skim")
    ap.add_argument("--cmssw_on_depot",
                    default="/depot/cms/users/shar1172/HZZ2l2nu_skim/CMSSW_14_0_2")
    ap.add_argument("--depot_bind", default="/depot/cms/users/shar1172",
                    help="writable /depot subtree to bind into the cmssw-el9 container")
    ap.add_argument("--redirector", default="root://xcache.cms.rcac.purdue.edu/")
    ap.add_argument("--account", default="cms")
    ap.add_argument("--partition", default="hammer-nodes")
    ap.add_argument("--time", default="04:00:00")
    ap.add_argument("--mem", default="4G")
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--max_parallel", type=int, default=200,
                    help="array throttle (%N)")
    ap.add_argument("--entries", type=int, default=0, help="post_proc --entriesToRun")
    ap.add_argument("--channels", default="2l2v", choices=["all", "4l", "2l2q", "2l2v"])
    ap.add_argument("--with_syst", action="store_true")
    ap.add_argument("--max_files", type=int, default=0,
                    help="debug: cap total input files")
    ap.add_argument("--dry_run", action="store_true", help="write files, do not sbatch")
    args = ap.parse_args()

    # ---- resolve inputs -------------------------------------------------------
    in_list = args.input_file
    if not os.path.isfile(in_list):
        in_list = os.path.join(PKG_DIR, "input_data_files", args.input_file)
    if not os.path.isfile(in_list):
        sys.exit("input list not found: %s" % args.input_file)

    proxy = find_proxy()
    if not proxy:
        sys.exit("no grid proxy found -- run `voms-proxy-init -voms cms`")

    if not os.path.isdir(args.cmssw_on_depot):
        sys.exit("CMSSW not staged on /depot (%s).\n  run: bash scripts/slurm/stage_to_depot.sh"
                 % args.cmssw_on_depot)

    with open(in_list) as fh:
        datasets = [l.strip() for l in fh
                    if l.strip() and not l.lstrip().startswith("#")]
    if not datasets:
        sys.exit("no datasets in %s" % in_list)

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    sub = "%s_%s" % (args.submission_name, ts)
    subdir = os.path.join(args.output_base, "submits", sub)
    skimbase = os.path.join(args.output_base, "skims", args.submission_name)
    os.makedirs(os.path.join(subdir, "logs"), exist_ok=True)

    # stash a copy of the proxy Slurm can read (/depot is visible, /tmp & ~ are not)
    depot_proxy = os.path.join(subdir, "x509_proxy")
    shutil.copy(proxy, depot_proxy)
    os.chmod(depot_proxy, 0o600)

    # ---- expand datasets ----------------------------------------------------
    tasks = []  # (idx, lfn_with_redirector, sample, outdir)
    for ds in datasets:
        tag = sample_tag(ds)
        files = das_files(ds)
        if not files:
            print("  !! no files: %s" % ds)
            continue
        outdir = os.path.join(skimbase, tag)
        os.makedirs(outdir, exist_ok=True)
        for f in files:
            # dasgoclient returns "/store/..."; xrootd wants  root://host//store/...
            lfn = args.redirector.rstrip("/") + "//" + f.lstrip("/")
            tasks.append((len(tasks), lfn, tag, outdir))
            if args.max_files and len(tasks) >= args.max_files:
                break
        print("  %-70s  %4d files" % (tag, len(files)))
        if args.max_files and len(tasks) >= args.max_files:
            break

    if not tasks:
        sys.exit("no input files resolved")

    tasks_tsv = os.path.join(subdir, "tasks.tsv")
    with open(tasks_tsv, "w") as fh:
        for t in tasks:
            fh.write("%d\t%s\t%s\t%s\n" % t)

    # ---- git provenance ---------------------------------------------------
    try:
        # infoCreaterGit reads os.environ['CMSSW_BASE'] unconditionally
        os.environ.setdefault("CMSSW_BASE", PKG_DIR.split("/src/")[0])
        sys.path.insert(0, os.path.join(PKG_DIR, "scripts", "utils"))
        import infoCreaterGit
        msg = "slurm submission %s" % sub
        info = infoCreaterGit.BasicInfoCreater(os.path.join(subdir, "summary.dat"), msg)
        info.generate_git_patch_and_log()
    except Exception as e:  # provenance is best-effort
        print("  (git provenance skipped: %s)" % e)

    # ---- the array job script ------------------------------------------
    njobs = len(tasks)
    job_sh = os.path.join(subdir, "job.sh")
    with open(job_sh, "w") as fh:
        fh.write("""#!/bin/bash
#SBATCH --job-name=hzzskim_{name}
#SBATCH --account={acct}
#SBATCH --partition={part}
#SBATCH --time={time}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --array=0-{last}%{par}
#SBATCH --output={sub}/logs/task_%a.out
#SBATCH --error={sub}/logs/task_%a.out

export HZZ_CMSSW="{cmssw}"
export HZZ_TASKS="{tasks}"
export HZZ_PROXY="{proxy}"
export HZZ_BIND="{bind}"
export HZZ_ENTRIES="{entries}"
export HZZ_CHANNELS="{channels}"
export HZZ_SYST="{syst}"

exec bash "{pkg}/scripts/slurm/job_wrapper.sh"
""".format(name=args.submission_name, acct=args.account, part=args.partition,
           time=args.time, cpus=args.cpus, mem=args.mem, last=njobs - 1,
           par=args.max_parallel, sub=subdir, cmssw=args.cmssw_on_depot,
           tasks=tasks_tsv, proxy=depot_proxy, bind=args.depot_bind,
           entries=args.entries, channels=args.channels,
           syst=("1" if args.with_syst else ""),
           pkg=os.path.join(args.cmssw_on_depot, PKG_REL)))
    os.chmod(job_sh, 0o755)

    print("\n%d array tasks  ->  %s" % (njobs, job_sh))
    print("skims  ->  %s/<sample>/" % skimbase)
    print("logs   ->  %s/logs/" % subdir)
    if args.dry_run:
        print("\n[dry-run] submit with:\n  sbatch %s" % job_sh)
        return
    r = subprocess.run(["sbatch", job_sh])
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
