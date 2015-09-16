import FWCore.ParameterSet.Config as cms

jetproperties = cms.EDProducer('JetProperties',
JetTag_               = cms.InputTag('slimmedJets'),
MinPt = cms.double(15),
doJEC = cms.bool(False),
RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
L1File = cms.string("PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt"),
L2File = cms.string("PHYS14_25_V2_All_L2Relative_AK4PFchs.txt"),
L3File = cms.string("PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt"),
#jecPayloadNames      = cms.vstring("JEC/PHYS14_25_V2::All_L3Absolute_AK4PFchs.txt","JEC/PHYS14_25_V2::All_L2Relative_AK4PFchs.txt","JEC/PHYS14_25_V2::All_L1FastJet_AK4PFchs.txt"),
#BTagInputTag	        = cms.string('combinedSecondaryVertexBJetTags'),
)
