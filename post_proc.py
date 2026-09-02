#!/usr/bin/env python3
import os
import sys
import argparse
import glob
import tempfile

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.met_phi_correction import METPhiCorrector, Campaign
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import (
    muonScaleRes2016pre,
    muonScaleRes2016,
    muonScaleRes2017,
    muonScaleRes2018,
    muonScaleRes2022,
    muonScaleRes2022EE,
    muonScaleRes2023,
    muonScaleRes2023BPix
)
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import createJMECorrector
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSFProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *

# Custom module imports
from modules.H4LCppModule import *
from modules.JetSFMaker import *
from modules.GenVarsProducer import *
from modules.keep_and_drop_list import keep_drop_rules_Data_MC, keep_drop_rules_GEN

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFile", default="", type=str, help="Input file name")
    parser.add_argument('-o', '--outputFile', default="skimmed_nano.root", type=str, help="Output file name")
    parser.add_argument('-outDir', '--outputDir', default=".", type=str, help="Output directory")
    parser.add_argument('-c', '--cutFlowFile', default="cutFlow.json", type=str, help="Cut flow file name")
    parser.add_argument("-n", "--entriesToRun", default=0, type=int, help="Set  to 0 if need to run over all entries else put number of entries to run")
    parser.add_argument("-d", "--DownloadFileToLocalThenRun", default=True, type=bool, help="Download file to local then run")
    parser.add_argument("--WithSyst", default=False, action="store_true", help="Do not run systematics")
    parser.add_argument("--DEBUG", default=False, action="store_true", help="Print debug information")
    parser.add_argument("--channels",  choices=["all", "4l", "2l2q", "2l2v"],  default="2l2v",
                        help="Channels to run: all, 4l, 2l2q, or 2l2v")
    return parser.parse_args()

def getListFromFile(filename):
    """Read file list from a text file."""
    with open(filename, "r") as file:
        return ["root://cms-xrd-global.cern.ch/" + line.strip() for line in file]

def create_temp_keep_drop_file(rules):
    """Create a temporary keep and drop file from a list of rules."""
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt')
    temp_file.write("\n".join(rules))
    temp_file.close()
    return temp_file.name

def main():
    args = parse_arguments()

    # Initial setup
    testfilelist = []
    modulesToRun = []
    isMC = True
    isFSR = True # set false for now
    isFiducialAna = True
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

    """Determine the year and type (MC or Data) of input ROOT file:
    For data the string "/data/" is always there. So, we take this
    as handle to decide if the root file is MC or data.
    """
    first_file = testfilelist[0]
    isMC = "/data/" not in first_file

    if "Summer22" in first_file or "Run2022" in first_file:
        """Summer22 and Run2022 for identification of 2022 MC and data respectiverly
        """
        year = 2022
        cfgFile = "config/Input_2022.yml"
        jsonFileName = "data/golden_json/Cert_Collisions2022_355100_362760_Golden.json"
        sfFileName = "DeepCSV_102XSF_V2.csv" # FIXME: Update for year 2022
        modulesToRun.extend([muonScaleRes2022()]) # FIXME: Update for year 2022
    if "UL18NanoAODv9" in first_file or "UL2018_MiniAODv2_NanoAODv9" in first_file:
        """UL2018 for identification of 2018 UL data and UL18 for identification of 2018 UL MC
        """
        year = 2018
        cfgFile = "config/Input_2018.yml"
        jsonFileName = "data/golden_json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2018()])
    if "UL17NanoAODv15" in first_file or "UL2017" in first_file:
        year = 2017
        cfgFile = "config/Input_2017.yml"
        jsonFileName="data/golden_json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2017()])

    if "20UL16NanoAODAPVv9" in first_file:
        year = 2016
        cfgFile = "config/Input_2016.yml"
        jsonFileName = "data/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2016pre()])
    if "20UL16NanoAODv9" in first_file:
        year = 2016
        cfgFile = "config/Input_2016.yml"
        jsonFileName = "data/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2016()])

    if "UL2018_NanoAODv15" in first_file or "UL18NanoAODv15" in first_file:
        year = 2018
        cfgFile = "config/Input_2018.yml"
        jsonFileName = "data/golden_json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2018()])


    #if cfgFile is None:
        #year = 2018   
        #cfgFile = "config/Input_2018.yml"
        #modulesToRun.extend([muonScaleRes2018()])

    H4LCppModule = lambda: HZZAnalysisCppProducer(year=year, cfgFile=cfgFile,
                                                  isMC=isMC, isFSR=isFSR,
                                                  cutFlowJSONFile=args.cutFlowFile,
                                                  channels=args.channels,
                                                  DEBUG=args.DEBUG
                                                  )
    print("systematic info: {}".format(args.WithSyst))
    print("Input json file: {}".format(jsonFileName))
    print("Input cfg file: {}".format(cfgFile))
    print("isMC: {}".format(isMC))
    print("isFSR: {}".format(isFSR))

    if isMC:
        # PU reweighting must run BEFORE H4LCppModule so that puWeight / puWeightUp /
        # puWeightDown exist when H4LCppModule builds overallEventWeight.
        # Use the UltraLegacy payloads (puWeight_UL20XX): they reweight to the UL
        # data pileup profile with a *fixed* MC profile (mcPileupUL20XX.root).
        # Do NOT use puAutoWeight_20XX here - those target the pre-UL (ReReco) data
        # profile and rebuild the MC profile from each input file, so the weight
        # depends on how files are split across jobs (non-reproducible).
        # NOTE: puWeight_UL2016 has no preVFP/postVFP split (nanoAOD-tools ships one
        # UL2016 file); a per-APV split needs the LUM correctionlib puWeights.json.gz.
        # FIXME: No PU weight for 2022 (needs the LUM Run-3 puWeights.json.gz) -
        # H4LCppModule then warns once and leaves puWeight out of overallEventWeight.
        if year == 2018: modulesToRun.extend([puWeight_UL2018()])
        if year == 2017: modulesToRun.extend([puWeight_UL2017()])
        if year == 2016: modulesToRun.extend([puWeight_UL2016()])

        GenVarModule = lambda : GenVarsProducer() # FIXME: Gen variable producer module is not working
        modulesToRun.extend([H4LCppModule(), GenVarModule()])
        #modulesToRun.extend([H4LCppModule()])
        if (args.WithSyst):
            # WARNING: this path only *stores* JES/JER-shifted jet/MET branches; it does
            # NOT re-run the C++ selection per variation, and createJMECorrector shifts
            # PF Type-I MET, not the PuppiMET the analysis uses. So JES/JER/unclustered-MET
            # shape systematics are NOT propagated to M_T end-to-end yet (needs an
            # H4LTools per-variation selection loop + PuppiMET shift branches).
            # jetType must be AK4PFchs for NanoAODv9 and AK4PFPuppi for the v15 re-nano.
            _ak4type = "AK4PFPuppi" if "NanoAODv15" in first_file else "AK4PFchs"
            jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = _ak4type)
            fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK8PFPuppi")
            # btag SF: the skim b-tags on DeepJet (Jet_btagDeepFlavB), so algo must be
            # "deepjet"; kept off until the b-veto in ZZSelection_2l2nu() is re-enabled.
            btagSF = lambda: btagSFProducer(era = "UL"+str(year), algo = "deepjet")
            # PU-jet-ID SF: Run 2 CHS (v9) only - reads Jet_puId, which the v15 re-nano
            # does not store. Do not enable for v15.
            puidSF = lambda: JetSFMaker("%s" % year)
            modulesToRun.extend([jetmetCorrector(), fatJetCorrector()])
            #modulesToRun.extend([jetmetCorrector(), fatJetCorrector(), puidSF()])   # v9 only
            #modulesToRun.extend([jetmetCorrector(), fatJetCorrector(), btagSF(), puidSF()])

        # INFO: Keep the `fwkJobReport=False` to trigger `haddnano.py`
        #            otherwise the output file will have larger size then expected. Reference: https://github.com/cms-nanoAOD/nanoAOD-tools/issues/249
        temp_keep_drop_file = create_temp_keep_drop_file(keep_drop_rules_GEN + keep_drop_rules_Data_MC)
        print("DEBUG: Keep and drop file: {}".format(temp_keep_drop_file))
        p=PostProcessor(args.outputDir,testfilelist, None, None,modules = modulesToRun,
                        provenance=True,fwkJobReport=True,
                        haddFileName=args.outputFile,
                        maxEntries=entriesToRun,
                        prefetch=DownloadFileToLocalThenRun, longTermCache= True,   # prefetch: download file to local then run, longTermCache: keep the file in local after running so that if it is present use local instead of downloading again
                        outputbranchsel=temp_keep_drop_file)
    else:
        modulesToRun.extend([H4LCppModule()])
        if (args.WithSyst):
            _ak4type = "AK4PFPuppi" if "NanoAODv15" in first_file else "AK4PFchs"
            jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = _ak4type)
            fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK8PFPuppi")
            modulesToRun.extend([jetmetCorrector(), fatJetCorrector()])
            #modulesToRun.extend([jetmetCorrector()])

        temp_keep_drop_file = create_temp_keep_drop_file(keep_drop_rules_Data_MC)
        print("DEBUG: Keep and drop file: {}".format(temp_keep_drop_file))
        p=PostProcessor(args.outputDir,testfilelist, None, None, modules = modulesToRun,
                        provenance=True, fwkJobReport=True,
                        haddFileName=args.outputFile,
                        jsonInput=jsonFileName,
                        maxEntries=entriesToRun,
                        prefetch=DownloadFileToLocalThenRun,  longTermCache= True,   # prefetch: download file to local then run, longTermCache: keep the file in local after running so that if it is present use local instead of downloading again
                        outputbranchsel=temp_keep_drop_file)

    p.run()

if __name__ == "__main__":
    main()
