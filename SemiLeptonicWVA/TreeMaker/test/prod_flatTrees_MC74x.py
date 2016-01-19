from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'WVGamma_MC_25ns'
config.section_('JobType')
config.JobType.psetName = 'TreeMaker/test/runMakeTreeFromMiniAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['Summer15_25nsV2_MC_L1FastJet_AK8PFchs.txt','Summer15_25nsV2_MC_L2Relative_AK8PFchs.txt','Summer15_25nsV2_MC_L3Absolute_AK8PFchs.txt','Summer15_25nsV2_MC_L1FastJet_AK4PFchs.txt','Summer15_25nsV2_MC_L2Relative_AK4PFchs.txt','Summer15_25nsV2_MC_L3Absolute_AK4PFchs.txt' ]
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.unitsPerJob = 1
config.Data.inputDBS = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'FileBased'
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
        
    config.General.requestName = 'WGamma_MC_25ns'
    config.JobType.pyCfgParams = ['global_tag=MCRUN2_74_V9::All', 'MC=True', 'name=RSGraviton600']
    config.Data.inputDataset = '/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
    config.Data.outLFNDirBase = '/store/user/jfaulkne/WG/'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()


    config.General.requestName = 'RSGraviton800'
    config.JobType.pyCfgParams = ['global_tag=MCRUN2_74_V9::All', 'MC=True', 'name=RSGraviton800']
    config.Data.inputDataset = '/RSGravToWWToLNQQ_kMpl01_M-800_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
    config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/RSGraviton800/'
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

