#!/usr/bin/env python
import os
import glob
import ROOT
import subprocess
import sys
import logging

logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')

def system(command):
    logging.debug("Executing command#12: {command}".format(command=command))
    #return subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)

def system_with_terminal_display(command):
    logging.info("Executing command#16: {command}".format(command=command))
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
        logging.warning("Zombie found: {}".format(fname))
    return isValid

def checkfaulty(fname, ref = None):
    if not ref:
        ref = ROOT.TFile.Open(fname)
    faultyfiles = []  # Initialize an empty list to hold the names of faulty files
    probe = ROOT.TFile.Open(fname)

    if not probe:
        logging.error("Could not open file {fname}".format(fname=fname))
        return False

    for e in ref.GetListOfKeys():
        name = e.GetName()
        try:
            otherObj = probe.GetListOfKeys().FindObject(name).ReadObj()
        except:
            faultyfiles.append(probe.GetName())

    probe.Close()

    if faultyfiles:
        logging.warning("Faulty files found: {faultyfiles}".format(faultyfiles=faultyfiles))
        return False

    return True

def isValidAndFaultFree(fname,ref = None):
    # First check if it's a valid ROOT file
    if not isValidRootFile(fname):
        return False

    # Then check for faulty keys
    return checkfaulty(fname, ref)

def merge_files(targetFile, filesToMerge):
    # If there are more than 100 files, split them into groups of 100
    if len(filesToMerge) > 100:
        logging.info('A lot of files to merge; this might take some time...')
        tempTargets = []
        tempFilesToMerge = [filesToMerge[x:x+100] for x in range(0, len(filesToMerge), 100)]

        temp_directory = "/tmp/rasharma/"
        for i, batch in enumerate(tempFilesToMerge):
            #tempTargetFile = targetFile.replace('.root', '-temp{}.root'.format(i))
            tempTargetFile = os.path.join(temp_directory, targetFile.replace('.root', '-temp{}.root'.format(i)))
            tempTargets.append(tempTargetFile)
            logging.info("tempTargetFile: {tempTargetFile}".format(tempTargetFile=tempTargetFile))
            logging.info("batch: {batch}".format(batch=batch))

            # Check if temporary target file already exists and is valid
            # if os.path.exists(tempTargetFile):
            #     if isValidAndFaultFree(tempTargetFile):
            #         continue
            #     else:
            #         logging.info("Removing temp hadd file {tempTargetFile}".format(tempTargetFile=tempTargetFile))
            #         system_with_terminal_display('rm {tempTargetFile}'.format(tempTargetFile=tempTargetFile))

            system_with_terminal_display('haddnano.py {0} {1}'.format(tempTargetFile, ' '.join(batch)))


        # Merge temporary files into the final target file
        targetFile = os.path.join(temp_directory, targetFile)
        system_with_terminal_display('haddnano.py {0} {1}'.format(targetFile, " ".join(tempTargets)))

        # Remove temporary files
        for tempTarget in tempTargets:
            system_with_terminal_display('rm {tempTarget}'.format(tempTarget=tempTarget))

    # If there are 100 or fewer files, directly merge them
    else:
        logging.info(system('haddnano.py {targetFile} {" ".join(filesToMerge)}'.format(targetFile=targetFile)))


def main():
    # Validate arguments
    if len(sys.argv) < 2:
        logging.error("Usage: script_name.py <submitVersion>")
        sys.exit(1)

    submitVersion = str(sys.argv[1])
    mainOutputDir = '/eos/user/a/avijay/HZZCondorjobRun2018_v9/DYJetsToLL_LHEFilterPtZ-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/231201_214931/{submitVersion}'.format(submitVersion=submitVersion)
    logging.info("submitVersion: {submitVersion}".format(submitVersion=submitVersion))
    logging.info("mainOutputDir: {mainOutputDir}".format(mainOutputDir=mainOutputDir))

    # List of already merged eras
    Alreadymerged = []

    eraDir = mainOutputDir
    # Loop through each era directory
    #for eraDir in glob.glob("{mainOutputDir}/*".format(mainOutputDir=mainOutputDir)):
        #logging.info("checking for extra folder")
        #if not os.path.isdir(eraDir):
    #if not os.path.exists(eraDir):
        #logging.info("not found extra folder")
        #continue
    era = os.path.basename(eraDir)
    logging.info("Processing era: {era}".format(era=era))

    #if era in Alreadymerged:
        #logging.info("Skipping this era as it is already merged.")
        #continue

    targetFile = os.path.join(eraDir, '.root').replace('/.root', '.root')
    targetFile = os.path.join(mainOutputDir.split('/')[-3], '.root').replace('/.root', '.root')
    filesToMerge = glob.glob(os.path.join(eraDir, '*.root'))
        #################
    logging.info("Before checking filesToMerge: Length of filesToMerge: {0}".format(len(filesToMerge)))
    logging.info("filesToMerge: {0}".format(filesToMerge))
        #################
    #if len(filesToMerge) == 0:
        #logging.info("No files to merge. Skipping.")
        #continue

        # Check if the target file already exists and is valid
    #if os.path.exists(targetFile):
        #if isValidAndFaultFree(targetFile):
            #logging.info("Seems hadd is already performed. Skipping.")
            #continue
        #else:
            #logging.info("Removing invalid target file {targetFile}".format(targetFile=targetFile))
            #system_with_terminal_display('rm {targetFile}'.format(targetFile=targetFile))

        # Merge the files
    logging.info("Before zombie: Length of filesToMerge: {0}".format(len(filesToMerge)))

    # filesToMerge = [f for f in filesToMerge if isValidAndFaultFree(f, ROOT.TFile.Open(filesToMerge[0]))]
    logging.info("After zombie: Length of filesToMerge: {0}".format(len(filesToMerge)))
    merge_files(targetFile, filesToMerge)

if __name__ == "__main__":
    main()
