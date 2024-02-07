#!/usr/bin/env python
import os
import sys
import argparse
import glob

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import createJMECorrector
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSFProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *

# Custom module imports
from H4Lmodule import *
from H4LCppModule import *
from JetSFMaker import *
# from GenVarsProducer import *

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFile", default="", type=str, help="Input file name")
    parser.add_argument("-n", "--entriesToRun", default=100, type=int, help="Set  to 0 if need to run over all entries else put number of entries to run")
    parser.add_argument("-d", "--DownloadFileToLocalThenRun", default=True, type=bool, help="Download file to local then run")
    parser.add_argument("--NOsyst", default=False, action="store_true", help="Do not run systematics")
    parser.add_argument("--DEBUG", default=False, action="store_true", help="Print debug information")
    return parser.parse_args()

def getListFromFile(filename):
    """Read file list from a text file."""
    with open(filename, "r") as file:
        return ["root://cms-xrd-global.cern.ch/" + line.strip() for line in file]

def main():
    args = parse_arguments()

    # Initial setup
    testfilelist = []
    modulesToRun = []
    isMC = True
    isFSR = True
    year = None
    cfgFile = None
    jsonFileName = None
    sfFileName = None

    entriesToRun = int(args.entriesToRun)
    DownloadFileToLocalThenRun = args.DownloadFileToLocalThenRun

    # Determine list of files to process
    if args.inputFile.endswith(".txt"):
        testfilelist = getListFromFile(args.inputFile)
    elif args.inputFile.endswith(".root"):
        testfilelist.append(args.inputFile)
    else:
        print("INFO: No input file specified. Using default file list.")
        testfilelist = getListFromFile("ExampleInputFileList.txt")
    print("DEBUG: Input file list: {}".format(testfilelist))
    if len(testfilelist) == 0:
        print("ERROR: No input files found. Exiting.")
        exit(1)

    # Determine the year and type (MC or Data)
    first_file = testfilelist[0]
    isMC = "/data/" not in first_file

    if "UL18" in first_file or "UL2018" in first_file:
        """UL2018 for identification of 2018 UL data and UL18 for identification of 2018 UL MC
        """
        year = 2018
        cfgFile = "Input_2018.yml"
        jsonFileName = "golden_Json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2018()])
    if "UL17" in first_file or "UL2017" in first_file:
        year = 2017
        cfgFile = "Input_2017.yml"
        jsonFileName="golden_Json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2017()])
    if "UL16" in first_file or "UL2016" in first_file:
        year = 2016
        jsonFileName = "golden_Json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2016()])
    H4LCppModule = lambda: HZZAnalysisCppProducer(year,cfgFile, isMC, isFSR, args.DEBUG)
    # GenVarModule = lambda : GenVarsProducer() # FIXME: Gen variable producer module is not working
    modulesToRun.extend([H4LCppModule()])
    # modulesToRun.extend([H4LCppModule(), GenVarModule()])
    # modulesToRun.extend([ GenVarModule()])

    print("systematic info: {}".format(args.NOsyst))
    print("Input json file: {}".format(jsonFileName))
    print("Input cfg file: {}".format(cfgFile))
    print("isMC: {}".format(isMC))
    print("isFSR: {}".format(isFSR))

    if isMC:
        if (not args.NOsyst):
            jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK4PFchs")
            fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK8PFPuppi")
            # btagSF = lambda: btagSFProducer("UL"+str(year), algo="deepjet",selectedWPs=['L','M','T','shape_corr'], sfFileName=sfFileName)
            btagSF = lambda: btagSFProducer(era = "UL"+str(year), algo = "deepcsv")
            puidSF = lambda: JetSFMaker("%s" % year)
            modulesToRun.extend([jetmetCorrector(), fatJetCorrector(), puidSF()])
        # # modulesToRun.extend([jetmetCorrector(), fatJetCorrector(), btagSF(), puidSF()])

        if year == 2018: modulesToRun.extend([puAutoWeight_2018()])
        if year == 2017: modulesToRun.extend([puAutoWeight_2017()])
        if year == 2016: modulesToRun.extend([puAutoWeight_2016()])

        p=PostProcessor(".",testfilelist, None, None,modules = modulesToRun,
                        provenance=True,fwkJobReport=False,
                        haddFileName="skimmed_nano_mc.root",
                        maxEntries=entriesToRun,
                        prefetch=DownloadFileToLocalThenRun, longTermCache= True,   # prefetch: download file to local then run, longTermCache: keep the file in local after running so that if it is present use local instead of downloading again
                        outputbranchsel="keep_and_drop.txt")
    else:
        if (not args.NOsyst):
            jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK4PFchs")
            fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK8PFPuppi")
            modulesToRun.extend([jetmetCorrector(), fatJetCorrector()])

        p=PostProcessor(".",testfilelist, None, None, modules = modulesToRun,
                        provenance=True, fwkJobReport=False,
                        haddFileName="skimmed_nano_data.root",
                        jsonInput=jsonFileName,
                        maxEntries=entriesToRun, prefetch=DownloadFileToLocalThenRun, longTermCache= True,
                        outputbranchsel="keep_and_drop_data.txt")

    p.run()

if __name__ == "__main__":
    main()

