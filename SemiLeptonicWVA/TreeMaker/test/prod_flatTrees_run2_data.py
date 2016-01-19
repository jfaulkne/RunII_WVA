from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'data_mu'
config.General.workArea = 'data_mu'
config.section_('JobType')
config.JobType.psetName = 'TreeMaker/test/runMakeTreeFromMiniAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['global_tag=74X_dataRun2_Prompt_v1', 'MC=False', 'isCrab=True']
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt','PHYS14_25_V2_All_L2Relative_AK8PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK8PFchs.txt','PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt','PHYS14_25_V2_All_L2Relative_AK4PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt' ]
config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2015B-PromptReco-v1/MINIAOD'
config.Data.unitsPerJob = 50
config.Data.lumiMask = 'json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
config.Data.inputDBS = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'LumiBased'
config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple_data_mu/'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand

    #Make sure you set this parameter (here or above in the config it does not matter)
    config.General.workArea = 'ntuple'

    def submit(config):
        res = crabCommand('submit', config = config)

    #########    From now on that's what users should modify: this is the a-la-CRAB2 configuration part.
        
    config.General.requestName = 'data_mu_prompt'
    config.Data.inputDataset = '/SingleMuon/Run2015B-PromptReco-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/data_mu_prompt/'
#    config.Data.lumiMask = 'json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'data_mu_17jul' #needed for runs<251585
    config.Data.inputDataset = '/SingleMuon/Run2015B-17Jul2015-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/data_mu_17jul/'
#    config.Data.lumiMask = 'json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'data_el_prompt'
    config.Data.inputDataset = '/SingleElectron/Run2015B-PromptReco-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/data_el_prompt/'
#    config.Data.lumiMask = 'json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'data_el_17jul' #needed for runs<251585
    config.Data.inputDataset = '/SingleElectron/Run2015B-17Jul2015-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/data_el_17jul/'
#    config.Data.lumiMask = 'json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    #...
