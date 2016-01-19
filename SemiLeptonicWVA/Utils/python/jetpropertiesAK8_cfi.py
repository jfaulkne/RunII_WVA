import FWCore.ParameterSet.Config as cms

jetpropertiesAK8 = cms.EDProducer('JetPropertiesAK8',
JetTag_               = cms.InputTag('slimmedJetsAK8'),
prunedJetTag_               = cms.InputTag('slimmedJetsAK8'),
softdropJetTag_               = cms.InputTag('slimmedJetsAK8'),
MinPt = cms.double(15),
doJEC = cms.bool(False),
doReclusteringForPrunedAndSoftdrop = cms.bool(False),
RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
BTagInputTag	        = cms.string('combinedSecondaryVertexBJetTags'),
L1File = cms.string("Summer15_25nsV5_DATA_L1FastJet_AK8PFchs.txt"),
L2File = cms.string("Summer15_25nsV5_DATA_L2Relative_AK8PFchs.txt"),
L3File = cms.string("Summer15_25nsV5_DATA_L3Absolute_AK8PFchs.txt"),
L2L3File = cms.string("Summer15_25nsV5_DATA_L2L3Residual_AK8PFchs.txt")
)
