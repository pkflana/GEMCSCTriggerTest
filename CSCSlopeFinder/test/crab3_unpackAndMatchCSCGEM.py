#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()
###2018runA  314472-318876
#section general
config.General.requestName = '2023D_Run370293_23Aug2023'
config.General.workArea = '2023D_Run370293_23Aug2023'#working dir 
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_unpack_CSCGEM_matcher.py'
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 1440 # 1440min = 24hours
config.JobType.numCores = 1
config.JobType.allowUndistributedCMSSW = True
#config.JobType.generator
#config.JobType.pyCfgParams
#config.JobType.inputFiles


config.Data.inputDataset = '/Muon1/Run2023D-ZMu-PromptReco-v1/RAW-RECO'
config.Data.runRange = '370293'

#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/daebi/GEMCSCTrigger/2023D_Run370293_23Aug2023/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.ignoreGlobalBlacklist = True
