import FWCore.ParameterSet.Config as cms

geneventinfo = cms.EDProducer('GenEventInfo',
geneventToken  = cms.InputTag("generator"),
)
