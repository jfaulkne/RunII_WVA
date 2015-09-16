import FWCore.ParameterSet.Config as cms

genJet = cms.EDProducer('GenJet',
  GenJetCollTag  = cms.InputTag("slimmedGenJets"),
)
