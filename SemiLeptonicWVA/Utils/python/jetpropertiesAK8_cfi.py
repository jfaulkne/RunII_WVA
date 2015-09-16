import FWCore.ParameterSet.Config as cms

jetpropertiesAK8 = cms.EDProducer('JetPropertiesAK8',
JetTag_               = cms.InputTag('slimmedJetsAK8'),
MinPt = cms.double(15),
#puppiJetTag_               = cms.InputTag('ak8PFJetsPuppi'),
doJEC = cms.bool(False),
RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
BTagInputTag	        = cms.string('combinedSecondaryVertexBJetTags'),
L1File = cms.string("PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt"),
L2File = cms.string("PHYS14_25_V2_All_L2Relative_AK8PFchs.txt"),
L3File = cms.string("PHYS14_25_V2_All_L3Absolute_AK8PFchs.txt"),
#jecPayloadNames      = cms.vstring("JEC/PHYS14_25_V2::All_L3Absolute_AK8PFchs.txt","JEC/PHYS14_25_V2::All_L2Relative_AK8PFchs.txt","JEC/PHYS14_25_V2::All_L1FastJet_AK8PFchs.txt"),
)
