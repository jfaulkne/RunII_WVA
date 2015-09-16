import FWCore.ParameterSet.Config as cms

metdouble = cms.EDProducer('METDouble',
METTag  = cms.InputTag("slimmedMETs"),
corrMet = cms.bool(False),
doJEC = cms.bool(False),
JetTag_               = cms.InputTag('slimmedJets'),
RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
L1File = cms.string("PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt"),
L2File = cms.string("PHYS14_25_V2_All_L2Relative_AK4PFchs.txt"),
L3File = cms.string("PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt"),
#jecPayloadNames      = cms.vstring("JEC/PHYS14_25_V2::All_L3Absolute_AK4PFchs.txt","JEC/PHYS14_25_V2::All_L2Relative_AK4PFchs.txt","JEC/PHYS14_25_V2::All_L1FastJet_AK4PFchs.txt"),
MuTag = cms.InputTag("slimmedMuons"),
)
