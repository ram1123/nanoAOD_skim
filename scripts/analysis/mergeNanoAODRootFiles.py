#!/usr/bin/env python3
import os
import glob
import ROOT
import subprocess
import sys
import logging
import argparse

logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')

def system(command):
    logging.debug("Executing command: {}".format(command))
    # return subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)

def system_with_terminal_display(command):
    logging.info("Executing command: {}".format(command))
    return subprocess.call(command, shell=True)

def isValidRootFile(fname):
    if not os.path.exists(fname):
        return False
    f = ROOT.TFile.Open(fname)
    if not f:
        return False
    isValid = not (f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered) or f.GetListOfKeys().IsEmpty())
    f.Close()
    if not isValid:
        logging.warning("Zombie file found: {}".format(fname))
    return isValid

def checkfaulty(fname, ref = None):
    if not ref:
        ref = ROOT.TFile.Open(fname)
    faultyfiles = []  # Initialize an empty list to hold the names of faulty files
    probe = ROOT.TFile.Open(fname)

    if not probe:
        logging.error("Could not open file {}".format(fname))
        return False

    for e in ref.GetListOfKeys():
        name = e.GetName()
        try:
            otherObj = probe.GetListOfKeys().FindObject(name).ReadObj()
        except:
            faultyfiles.append(probe.GetName())

    probe.Close()

    if faultyfiles:
        logging.warning("Faulty files found: {}".format(", ".join(faultyfiles)))
        return False

    return True

def isValidAndFaultFree(fname, ref=None):
    # First check if it's a valid ROOT file
    if not isValidRootFile(fname):
        return False

    # Then check for faulty keys
    return checkfaulty(fname, ref)

def merge_files(targetFile, filesToMerge):
    logging.info("Merging {} files into: {}".format(len(filesToMerge), targetFile))
    if len(filesToMerge) > 100:
        logging.info('A lot of files to merge; this might take some time...')
        tempTargets = []
        tempFilesToMerge = [filesToMerge[x:x+100] for x in range(0, len(filesToMerge), 100)]

        temp_directory = os.path.join("/tmp", os.getenv('USER', 'default_user'))
        logging.debug("Using temporary directory: {}".format(temp_directory))
        for i, batch in enumerate(tempFilesToMerge):
            tempTargetFile = os.path.join(temp_directory, os.path.basename(targetFile).replace('.root', '-temp{}.root'.format(i)))
            logging.debug("Merging batch {0} into temp file {1}".format(i, tempTargetFile))
            tempTargets.append(tempTargetFile)
            # Check if temporary target file already exists and is valid
            if os.path.exists(tempTargetFile):
                if isValidAndFaultFree(tempTargetFile):
                    continue
                else:
                    logging.info("Removing temp hadd file {tempTargetFile}".format(tempTargetFile=tempTargetFile))
                    system_with_terminal_display('rm {tempTargetFile}'.format(tempTargetFile=tempTargetFile))

            system_with_terminal_display('haddnano.py {0} {1}'.format(tempTargetFile, ' '.join(batch)))

        # Final merge
        system_with_terminal_display('haddnano.py {0} {1}'.format(targetFile, " ".join(tempTargets)))
        # Cleanup
        for tempTarget in tempTargets:
            logging.debug("Removing temp hadd file {tempTarget}".format(tempTarget=tempTarget))
            os.remove(tempTarget)
    else:
        logging.info("Files are < 100; merging directly.")
        logging.info("haddnano.py {0} {1}".format(targetFile, " ".join(filesToMerge)))
        system_with_terminal_display('haddnano.py {0} {1}'.format(targetFile, " ".join(filesToMerge)))

def main():
    parser = argparse.ArgumentParser(description="Merge ROOT files using haddnano.py.")
    parser.add_argument('-i', '--inputDir', type=str, required=True, help="Path of the input directory that contains ROOT files to be merged.")
    parser.add_argument('-o', '--outputDir', type=str, required=True, help="Path of the output directory where the merged ROOT file will be saved.")
    parser.add_argument('-f', '--outputFile', type=str, required=True, help="Name of the hadd-ed output ROOT file.")
    parser.add_argument('-r', '--recursive', action='store_true', help="Search for ROOT files recursively in the input directory.")

    args = parser.parse_args()

    inputDir = args.inputDir
    outputDir = args.outputDir
    outputFile = args.outputFile

    logging.info("Input directory: {}".format(inputDir))
    logging.info("Output directory: {}".format(outputDir))
    logging.info("Output file: {}".format(outputFile))

    if not os.path.isdir(inputDir):
        logging.error("The specified input directory does not exist: {}".format(inputDir))
        sys.exit(1)

    targetFile = os.path.join(outputDir, outputFile)

    filesToMerge = []
    if args.recursive:
        logging.info("Searching for ROOT files recursively in the input directory.")
        # filesToMerge = glob.glob(os.path.join(inputDir, '**', '*.root'), recursive=True) # Works only with Python 3
        for root, dirnames, filenames in os.walk(inputDir):
            for filename in glob.glob(os.path.join(root, '*.root')):
                filesToMerge.append(os.path.join(root, filename))
    else:
        logging.info("Searching for ROOT files in the input directory only.")
        filesToMerge = glob.glob(os.path.join(inputDir, '*.root')) # Does not search recursively


    if not filesToMerge:
        logging.error("No ROOT files found in the specified directory.")
        sys.exit(1)

    logging.info("Found {} ROOT files to merge.".format(len(filesToMerge)))
    logging.info("Merging files into: {}".format(targetFile))
    logging.debug("Files to merge: {}".format(filesToMerge))

    # Check if the target file already exists and is valid
    if os.path.exists(targetFile):
        if isValidAndFaultFree(targetFile):
            logging.info("The target file already exists and is valid.")
            sys.exit(0)
        else:
            logging.info("Removing existing target file: {}".format(targetFile))
            system_with_terminal_display('rm {targetFile}'.format(targetFile=targetFile))

    # Check the validity of the files to merge
    validFiles = []
    for f in filesToMerge:
        if isValidAndFaultFree(f):
            validFiles.append(f)
    filesToMerge = validFiles

    merge_files(targetFile, filesToMerge)

if __name__ == "__main__":
    main()
