import FWCore.ParameterSet.Config as cms

jetproperties = cms.EDProducer('JetProperties',
JetTag_               = cms.InputTag('slimmedJets'),
MinPt = cms.double(15),
doJEC = cms.bool(False),
RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
L1File = cms.string("Summer15_25nsV5_DATA_L1FastJet_AK4PFchs.txt"),
L2File = cms.string("Summer15_25nsV5_DATA_L2Relative_AK4PFchs.txt"),
L3File = cms.string("Summer15_25nsV5_DATA_L3Absolute_AK4PFchs.txt"),
L2L3File = cms.string("Summer15_25nsV5_DATA_L2L3Residual_AK4PFchs.txt")
)
