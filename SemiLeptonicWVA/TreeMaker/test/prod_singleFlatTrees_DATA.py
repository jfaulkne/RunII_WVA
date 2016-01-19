from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'data2'
config.General.workArea = 'data2'
config.section_('JobType')
config.JobType.psetName = 'TreeMaker/test/runMakeTreeFromMiniAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['global_tag=GR_R_74_V9A::All', 'MC=False']
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt','PHYS14_25_V2_All_L2Relative_AK8PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK8PFchs.txt','PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt','PHYS14_25_V2_All_L2Relative_AK4PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt' ]
config.section_('Data')
config.Data.inputDataset = '/Commissioning/lbrianza-crab_data_Commissioning_Run2015A-PromptReco-v1-e050d919924aa3d4cdb5e91105bf430f/USER'
config.Data.unitsPerJob = 50
config.Data.inputDBS = 'phys03' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'LumiBased'
config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple_data2/'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
