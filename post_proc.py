#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from H4Lmodule import *
from H4LCppModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *

ifRunningOnCondor = False

testfilelist = []

if ifRunningOnCondor:
    testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/120000/3ED05633-EBB7-4A44-8F9D-CD956490BCFD.root"
    testfilelist.append(testfile)

else:
    if len(sys.argv) > 1:
       InputFileList = sys.argv[1]
    else:
       InputFileList = "ExampleInputFileList.txt"
    with open(InputFileList, 'r') as file:
      for line in file:
        # Remove newline characters
        line = line.strip()
        # Append the line to the list with the "root://cms-xrd-global.cern.ch//" prefix
        testfilelist.append("root://cms-xrd-global.cern.ch/" + line)

# Set  entriesToRun = 0 if need to run over all entries else put number of entries to run
entriesToRun = 0 if ifRunningOnCondor else 100

isMC = True
isFSR = False

outputbranchsel = "keep_and_drop.txt"
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun=True

if testfilelist[0].find("/data/") != -1:
    isMC = False
    outputbranchsel = "keep_and_drop_data.txt"

if testfilelist[0].find("UL18") != -1 or testfilelist[0].find("UL2018") != -1: # UL2018 for identification of 2018 UL data and UL18 for identification of 2018 UL MC
    year = 2018
    cfgFile = 'Input_2018.yml'

if testfilelist[0].find("UL17") != -1 or testfilelist[0].find("UL2017") != -1:
    year = 2017
    cfgFile = 'Input_2017.yml'

H4LCppModule = lambda: HZZAnalysisCppProducer(year,cfgFile, isMC, isFSR)
print("outputbranchsel = ",outputbranchsel)
print("isMC = ",isMC)
print("isFSR = ",isFSR)
#p=PostProcessor(".",[testfile],"",None,[H4LCppModule()],provenance=True,fwkJobReport=False,haddFileName="nano_M125.root",maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun,outputbranchsel="keep_and_drop.txt")
p=PostProcessor(".",testfilelist,"",None,[H4LCppModule()],provenance=True,fwkJobReport=False,haddFileName="nano_M125_cpp.root",maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun,outputbranchsel=outputbranchsel)

p.run()
print "DONE"
