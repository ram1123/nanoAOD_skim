#!/usr/bin/env python
import os,sys

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from H4Lmodule import *
from H4LCppModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *

ifRunningOnCondor = False
WhichSample = "2l2q" # options: "4l", "2l2q", "2l2nu"

testfilelist = []

if ifRunningOnCondor:
    testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/120000/3ED05633-EBB7-4A44-8F9D-CD956490BCFD.root"
    testfilelist.append(testfile)

if WhichSample == "4l" and (not ifRunningOnCondor):
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/130000/4743B911-1EA3-7E46-959B-93F466ED622F.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/130000/DE3F93DE-DAF6-DD4A-AC3B-13ACE92EAD68.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/2430000/336F0076-347C-AC40-B1B4-E31E14323B81.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/260000/100D0E96-C700-5244-874F-38DB4B410228.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/260000/28512C40-9C25-864A-98E9-4746C0471E63.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/260000/2BAB6B06-FFE2-DC4F-A57F-3A33B3326A30.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/260000/453FC70D-8164-1B4B-A756-81E7432C1D61.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/260000/5C1E0BC2-F153-E946-9B61-470E8AF85C58.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/260000/9DB74D60-0434-A04C-B430-173D6D6538C3.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/260000/DB3B1648-246F-314E-9DD8-FACEA0AE62F6.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/260000/FD91C4F0-7DE3-4947-817D-3EA957A0BC50.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/270000/6021CDF7-AE6F-464D-81E2-AB75ABC4809D.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/270000/77B90DAA-35B9-004E-9C49-52E516FB15D6.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/270000/C01084F4-328D-0146-B71E-B167AB6A7E86.root")

# 2l2q ggH Sample 1000 GeV

#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
#     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/3822B676-3F8E-5F46-96E2-1BB4040641B9.root"
## BKG
     testfile = "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2540000/A58AE71A-0273-8546-8FBD-3F729CAF6516.root"
if WhichSample == "2l2q"  and (not ifRunningOnCondor):
    ## UL2018
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/130000/39E6146D-1A47-5042-B843-5F7263B262E9.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/270000/A226A60B-9BD4-B64B-B53E-A26DB011BD7B.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/270000/2E79070C-3860-C84E-A0E9-B63E85E42B8A.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/120000/0188520D-6137-9546-8D8A-546D089E4F79.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/280000/90E68AE8-1691-CC40-B817-71AF5A46B429.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/270000/F2D8B29B-E543-2444-86BB-726AF4060A2A.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/280000/7F817CD1-84C1-5347-AEEA-F069CAE84799.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/120000/CB8568FB-F94F-5C4D-B763-2390FE18E0C9.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/270000/F2261826-B5D4-A544-9AF0-03FCF04EDFCE.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/280000/844A1B6D-AC02-3641-B9D2-08446FBEA30B.root")
#    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/F920B631-39E7-EE4C-9D8D-B242A2103115.root")
    testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2540000/A58AE71A-0273-8546-8FBD-3F729CAF6516.root")
    ##testfilelist.append("root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9_GT36-v1/2820000/06EFD7A6-A6F1-E145-B1F8-C6794FE4EDDD.root")

    ## UL2017
    # testfilelist.append("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/40000/4E6C744B-023A-E14A-9ADD-663C7101E98E.root")

# 2l2nu ggH Sample 1000 GeV
if WhichSample == "2l2nu"  and (not ifRunningOnCondor): #2l2nu
    testfilelist.append("root://cms-xrd-global.cern.ch/")

entriesToRun = 100 # 0 if need to run over all entries else put number of entries to run
isMC = True
# Keep DownloadFileToLocalThenRun=True this should reduce the file read error from eos.
DownloadFileToLocalThenRun=True

if testfilelist[0].find("/data/") != -1: isMC = False
if testfilelist[0].find("UL18") != -1:
    year = 2018
    cfgFile = 'Input_2018.yml'

if testfilelist[0].find("UL17") != -1:
    year = 2017
    cfgFile = 'Input_2017.yml'
H4LCppModule = lambda: HZZAnalysisCppProducer(isMC, year,cfgFile)
#p=PostProcessor(".",[testfile],"",None,[H4LCppModule()],provenance=True,fwkJobReport=False,haddFileName="nano_M125.root",maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun,outputbranchsel="keep_and_drop.txt")
p=PostProcessor(".",testfilelist,"",None,[H4LCppModule()],provenance=True,fwkJobReport=False,haddFileName="nano_M125_cpp.root",maxEntries=entriesToRun,prefetch=DownloadFileToLocalThenRun,outputbranchsel="keep_and_drop.txt")

p.run()
print "DONE"
