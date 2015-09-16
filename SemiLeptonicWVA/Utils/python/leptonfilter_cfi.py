import FWCore.ParameterSet.Config as cms

leptonFilter = cms.EDFilter("LeptonFilter",
                           electronsInputTag = cms.InputTag("slimmedElectrons"),
                           muonsInputTag = cms.InputTag("slimmedMuons"),
                           eleFilterPtCut = cms.double(10.0),
                           muFilterPtCut = cms.double(5.0)
                           )
