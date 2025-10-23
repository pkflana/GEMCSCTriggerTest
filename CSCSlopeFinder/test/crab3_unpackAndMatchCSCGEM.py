#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()
###2018runA  314472-318876
#section general


config.General.requestName = '2024_ZMu_31Mar2025' #'2024G_ZeroBias_22Aug2024' #'2024E_Muon0_Zmu_28June2024'
config.General.workArea = '2024_ZMu_31Mar2025' #'2024G_ZeroBias_22Aug2024'#working dir

#config.General.requestName = '2024_ZeroBias_31Mar2025'
#config.General.workArea = '2024_ZeroBias_31Mar2025'#working dir

config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_unpack_CSCGEM_matcher.py'
config.JobType.maxMemoryMB = 4000
config.JobType.maxJobRuntimeMin = 420 # 7hrs # 1440min = 24hours
config.JobType.numCores = 1
config.JobType.allowUndistributedCMSSW = True
# config.JobType.inputFiles = ['/afs/cern.ch/work/d/daebi/gemcsctrigger/CMSSW_14_0_7/src/GEMCSCTriggerTest/CSCSlopeFinder/luts']
# config.JobType.inputFiles = ['/afs/cern.ch/work/v/vdamante/GEM_CSC_Development/GEMCSCTriggerTest/CSCSlopeFinder/luts']

config.JobType.inputFiles = ['/afs/cern.ch/work/v/vdamante/GEM_CSC_Development/GEMCSCTriggerTest/soft/CMSSW_14_2_0_pre1/src/GEMCSCTriggerTest/CSCSlopeFinder/luts']

#config.JobType.generator
#config.JobType.pyCfgParams
#config.JobType.inputFiles


# config.Data.inputDataset = '/EphemeralZeroBias0/Run2024I-PromptReco-v1/MINIAOD'
# config.Data.secondaryInputDataset = '/EphemeralZeroBias0/Run2024I-v1/RAW'
# config.Data.useParent = True # This allows us to put MiniAOD as the input, and it will find the parent files for the RAW parts
# # I never got the Ephemeral to work ):

# config.Data.inputDataset = '/ZeroBias/Run2025C-LogError-PromptReco-v1/RAW-RECO'
# config.Data.inputDataset = '/Muon0/Run2024I-ZMu-PromptReco-v1/RAW-RECO'
# config.Data.inputDataset = '/Muon0/Run2025C-ZMu-PromptReco-v1/RAW-RECO'

config.Data.inputDataset = '/SingleMuon_Pt-2To50_Eta-0p9To2p85-gun/Phase2Spring24DIGIRECOMiniAOD-PU140_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD'
#config.Data.runRange = '381384'

#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
# config.Data.splitting = 'LumiBased'
# config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 5
#config.Data.outLFNDirBase = '/store/user/daebi/GEMCSCTrigger/2024_ZeroBias/'
config.Data.outLFNDirBase = '/store/user/vdamante/GEMCSCTrigger/2024_ZMu/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.ignoreGlobalBlacklist = True
