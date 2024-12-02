from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = Configuration()

config.section_("General")
config.General.requestName = 'nanoAOD_testing_2'
config.General.workArea = '/uscms/home/rasharma/nobackup/nanoAOD/CMSSW_10_2_22/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/crab/crab_projects_sendPythonFolderFalse'
config.General.transferLogs=True
config.General.transferOutputs = True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['crab_script.py','../../../../../scripts/haddnano.py','keep_and_drop.txt'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder     = True
config.JobType.allowUndistributedCMSSW = True
config.section_("Data")
config.Data.inputDataset = '/WplusTo2JWminusTo2LJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM'
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/user/rasharma/NanoPostTemp'
#config.Data.outLFNDirBase = '/store/user/%s/NanoPostTemp' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.ignoreLocality = False
config.Data.outputDatasetTag = 'NanoTestPost'
config.section_("Site")
config.Site.storageSite = "T3_US_FNALLPC"
