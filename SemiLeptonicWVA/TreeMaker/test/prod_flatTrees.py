from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'RSGraviton1000'
config.General.workArea = 'RSGraviton1000'
config.section_('JobType')
config.JobType.psetName = 'TreeMaker/test/runMakeTreeFromMiniAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['global_tag=PHYS14_25_V2::All', 'MC=True']
config.JobType.inputFiles = ['PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt','PHYS14_25_V2_All_L2Relative_AK8PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK8PFchs.txt','PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt','PHYS14_25_V2_All_L2Relative_AK4PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt' ]
config.JobType.allowNonProductionCMSSW = True
config.section_('Data')
config.Data.inputDataset = '/RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v2/MINIAODSIM'
config.Data.unitsPerJob = 1
config.Data.inputDBS = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'FileBased'
config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/RSGraviton1000/'
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
    config.General.requestName = 'RSGraviton1000'
    config.Data.inputDataset = '/RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v2/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton1000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/RSGraviton1000/'
#    config.Data.unitsPerJob = 2
#    config.Data.totalUnits = 4
    submit(config)
    
    config.General.requestName = 'RSGraviton4000'
    config.Data.inputDataset = '/RSGravitonToWW_kMpl01_M_4000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/RSGraviton4000/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'WJets100'
    config.Data.inputDataset = '/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/WJets100/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'WJets200'
    config.Data.inputDataset = '/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/WJets200/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'WJets400'
    config.Data.inputDataset = '/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/WJets400/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'WJets600'
    config.Data.inputDataset = '/WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/WJets600/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'TTbar'
    config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/TTbar/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'sch'
    config.Data.inputDataset = '/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/sch/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'tch'
    config.Data.inputDataset = '/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/tch/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'tWch'
    config.Data.inputDataset = '/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/tWch/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'sch_bar'
    config.Data.inputDataset = '/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/sch_bar/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)
    
    
    config.General.requestName = 'tch_bar'
    config.Data.inputDataset = '/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/tch_bar/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)
 
    
    config.General.requestName = 'tWch_bar'
    config.Data.inputDataset = '/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#    config.General.workArea = 'ntuple/RSGraviton4000'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/tWch_bar/'
#    config.Data.unitsPerJob = 3
#    config.Data.totalUnits = 8
    submit(config)

    config.General.requestName = 'RSGraviton2000'
    config.Data.inputDataset = '/RSGravitonToWW_kMpl01_M_2000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/RSGraviton2000/'
    submit(config)

    config.General.requestName = 'RSGraviton3000'
    config.Data.inputDataset = '/RSGravitonToWW_kMpl01_M_3000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
    config.Data.outLFN = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple/RSGraviton3000/'
    submit(config)
    

    #...
