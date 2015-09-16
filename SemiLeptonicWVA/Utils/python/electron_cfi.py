import FWCore.ParameterSet.Config as cms

electron = cms.EDProducer('Electron',
VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
EleTag = cms.InputTag("slimmedElectrons"),
MinPt = cms.double(15),
RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
ConversionTag = cms.InputTag('reducedEgamma:reducedConversions'),
BeamSpotTag = cms.InputTag('offlineBeamSpot'),
eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60")
)
