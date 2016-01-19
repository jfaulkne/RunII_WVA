import FWCore.ParameterSet.Config as cms

genPhotonRecoCand = cms.EDProducer('GenPhotonRecoCand',
  PrunedGenParticleTag  = cms.InputTag("prunedGenParticles"),
)
