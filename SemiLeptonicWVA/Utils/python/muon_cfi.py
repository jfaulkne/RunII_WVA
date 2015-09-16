import FWCore.ParameterSet.Config as cms

muon = cms.EDProducer('Muon',
VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
MuTag = cms.InputTag("slimmedMuons"),
MinPt = cms.double(15),
RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
)
