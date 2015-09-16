import FWCore.ParameterSet.Config as cms

genJetAK8 = cms.EDProducer('GenJetAK8',
  GenJetCollTag  = cms.InputTag("selectedGenJetsAK8"),
)
