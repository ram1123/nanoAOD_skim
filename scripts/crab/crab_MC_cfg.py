from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = Configuration()

config.section_("General")
config.General.requestName = 'NanoPost1'
#config.General.workArea = 'crab_project'
config.General.workArea = '/uscmst1b_scratch/lpc1/3DayLifetime/rasharma/crab_project'
config.General.transferLogs=True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['crab_scriptMC.py','../../../../../scripts/haddnano.py','keep_and_drop.txt'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder     = True
config.JobType.allowUndistributedCMSSW = True
config.section_("Data")
config.Data.inputDataset = '/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM'
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 10

config.Data.outLFNDirBase = '/store/user/rasharma/NanoAOD_Skimed/2018'
#config.Data.outLFNDirBase = '/store/user/%s/NanoPostTemp' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'NanoTestPost'
config.section_("Site")
config.Site.storageSite = "T3_US_FNALLPC"
