from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'data_mu'
config.General.workArea = 'data_mu'
config.section_('JobType')
config.JobType.psetName = 'TreeMaker/test/runMakeTreeFromMiniAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['global_tag=74X_dataRun2_Prompt_v2', 'MC=False', 'isCrab=True', 'doJECCorrection=True']
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['Summer15_25nsV2_DATA_L1FastJet_AK8PFchs.txt','Summer15_25nsV2_DATA_L2Relative_AK8PFchs.txt','Summer15_25nsV2_DATA_L3Absolute_AK8PFchs.txt','Summer15_25nsV2_DATA_L1FastJet_AK4PFchs.txt','Summer15_25nsV2_DATA_L2Relative_AK4PFchs.txt','Summer15_25nsV2_DATA_L3Absolute_AK4PFchs.txt' ]
config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2015B-PromptReco-v1/MINIAOD'
config.Data.unitsPerJob = 50
config.Data.lumiMask = 'json/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
config.Data.inputDBS = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'LumiBased'
config.Data.outLFNDirBase = '/store/user/jfaulkn3/Data_Run2015CD/'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand

    #Make sure you set this parameter (here or above in the config it does not matter)
    config.General.workArea = 'ntuple'

    def submit(config):
        res = crabCommand('submit', config = config)

    #########    From now on that's what users should modify: this is the a-la-CRAB2 configuration part.
        
    config.General.requestName = 'data_mu_prompt_25ns_runC'
    config.Data.inputDataset = '/SingleMuon/Run2015C-PromptReco-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/jfaulkn3/Data_Run2015CD/data_mu_prompt_25ns_runC/'
#    config.Data.lumiMask = 'json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'data_el_prompt_25ns_runC'
    config.Data.inputDataset = '/SingleElectron/Run2015C-PromptReco-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/jfaulkn3/Data_Run2015CD/data_el_prompt_25ns_runC/'
#    config.Data.lumiMask = 'json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'data_mu_prompt_25ns_runD'
    config.Data.inputDataset = '/SingleMuon/Run2015D-PromptReco-v3/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/jfaulkn3/Data_Run2015CD/data_mu_prompt_25ns_runD/'
#    config.Data.lumiMask = 'json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'data_el_prompt_25ns_runD'
    config.Data.inputDataset = '/SingleElectron/Run2015D-PromptReco-v3/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/jfaulkn3/Data_Run2015CD/data_el_prompt_25ns_runD/'
#    config.Data.lumiMask = 'json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    #...
