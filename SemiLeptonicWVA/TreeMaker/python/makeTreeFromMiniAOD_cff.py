# $Id: makeTreeFromPAT_cff.py,v 1.16 2013/01/24 15:42:53 mschrode Exp $
#

import FWCore.ParameterSet.Config as cms
import os

def makeTreeTreeFromMiniAOD(process,
outFileName,
NJetsMin=2,
HTMin=350.,
MHTMin=0.,
reportEveryEvt=10,
testFileName="",
Global_Tag="",
METFiltersProcess="",
MC=False,
debug = False,
QCD=False,
LostLepton=False,
numProcessedEvt=1000,
doAK8Reclustering=False,
doJECCorrection=False,
doPuppi=False,
leptonFilter=True,
genJetsAK8Reclustering=True,
customizeHBHENoiseForEarlyData=False,
customizeHBHENoiseForRun2015D=True,
jsonFileName="",
reDoPruningAndSoftdrop=False,
isCrab=False):

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
    if (MC):
        customizeHBHENoiseForRun2015D=False

    process.GlobalTag.globaltag = Global_Tag

    ## Added Geometry cfi files ###
    ## Not in current EXO-WW configuration ##
    process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
    process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
    process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

    ## --- Log output ------------------------------------------------------
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.MessageLogger.cerr = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
        )
    process.MessageLogger.cout = cms.untracked.PSet(
        INFO = cms.untracked.PSet(reportEvery = cms.untracked.int32(reportEveryEvt))
        )
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        ) 


    ## --- Files to process ------------------------------------------------
    process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(numProcessedEvt)
        )
    process.source = cms.Source("PoolSource",

        fileNames = cms.untracked.vstring(testFileName)
        )

 ## ----------------------------------------------------------------------------------------------
## Triggers
## ----------------------------------------------------------------------------------------------
# The trigger results are saved to the tree as a vector
# Three vectors are saved:
# 1) names of the triggers
# 2) trigger results
# 3) trigger prescales
# the indexing of these vectors must match
# If the version number of the input trigger name is omitted,
# any matching trigger will be included (default behavior)
    from SemiLeptonicWVA.Utils.triggerproducer_cfi import triggerProducer
    process.TriggerProducer = triggerProducer.clone( 
        trigTagArg1 = cms.string('TriggerResults'),
        trigTagArg2 = cms.string(''),
        trigTagArg3 = cms.string('HLT'),
        prescaleTagArg1 = cms.string('patTrigger'),
        prescaleTagArg2 = cms.string(''),
        prescaleTagArg3 = cms.string(''),
        triggerNameList = cms.vstring( # list of trigger names
            'HLT_Ele22_eta2p1_WPTight_Gsf_v',
            'HLT_IsoMu20_v',
            'HLT_IsoTkMu20_v',
	    'HLT_Ele27_WP85_Gsf_v'
            )
        )

    ## --- Output file -----------------------------------------------------
    process.TFileService = cms.Service(
        "TFileService",
        fileName = cms.string(outFileName+".root")
        )

    ############### JSON Filter            
    import FWCore.PythonUtilities.LumiList as LumiList
    import sys

    if not MC:
        if(len(jsonFileName)>0):
            import FWCore.PythonUtilities.LumiList as LumiList
            process.source.lumisToProcess = LumiList.LumiList(filename = jsonFileName).getVLuminosityBlockRange()
        else:
            print "ERROR!! running on data with no json file applied!"
            sys.exit()

####### some gen infos
    from SemiLeptonicWVA.Utils.geneventinfo_cfi import geneventinfo
    process.GenEventInfo = geneventinfo.clone()
    
###Lepton Filter
    process.filterSeq = cms.Sequence ()

    process.load('SemiLeptonicWVA.Utils.leptonfilter_cfi')

    process.leptonFilter.electronsInputTag = cms.InputTag("slimmedElectrons")
    process.leptonFilter.muonsInputTag = cms.InputTag("slimmedMuons")
    process.leptonFilter.eleFilterPtCut = cms.double(20.0)
    process.leptonFilter.muFilterPtCut = cms.double(20.0)
    
    if (leptonFilter):
        process.filterSeq = cms.Sequence (process.leptonFilter)
    
       ## --- Setup of TreeMaker ----------------------------------------------
    FilterNames = cms.VInputTag()
    FilterNames.append(cms.InputTag("HBHENoiseFilterRA2","HBHENoiseFilterResult","PAT"))
    FilterNames.append(cms.InputTag("beamHaloFilter"))
    FilterNames.append(cms.InputTag("eeNoiseFilter"))
    FilterNames.append(cms.InputTag("trackingFailureFilter"))
    FilterNames.append(cms.InputTag("inconsistentMuons"))
    FilterNames.append(cms.InputTag("greedyMuons"))
    FilterNames.append(cms.InputTag("ra2EcalTPFilter"))
    FilterNames.append(cms.InputTag("ra2EcalBEFilter"))
    FilterNames.append(cms.InputTag("hcalLaserEventFilter"))
    FilterNames.append(cms.InputTag("ecalLaserCorrFilter"))
    FilterNames.append(cms.InputTag("eeBadScFilter"))
    FilterNames.append(cms.InputTag("PBNRFilter"))
    FilterNames.append(cms.InputTag("HCALLaserEvtFilterList2012"))
    FilterNames.append(cms.InputTag("manystripclus53X"))
    FilterNames.append(cms.InputTag("toomanystripclus53X"))
    FilterNames.append(cms.InputTag("logErrorTooManyClusters"))
    FilterNames.append(cms.InputTag("RA2HONoiseFilter"))

    
    ## --- Setup WeightProducer -------------------------------------------
    from SemiLeptonicWVA.WeightProducer.getWeightProducer_cff import getWeightProducer
    process.WeightProducer = getWeightProducer(testFileName)
    process.WeightProducer.Lumi                       = cms.double(5000)
    process.WeightProducer.PU                         = cms.int32(0) # PU S10 3 for S10 2 for S7
    process.WeightProducer.FileNamePUDataDistribution = cms.string("NONE")
    print process.WeightProducer.PU

    from RecoBTag.Configuration.RecoBTag_cff import *
    from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import *
    process.slimmedJetsPFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
      j2tParametersVX,
      jets = cms.InputTag("iterativeCone5PFJets")
    )
    process.slimmedJetsPFJetTracksAssociatorAtVertex.jets = "slimmedJets"
    process.slimmedJetsPFJetTracksAssociatorAtVertex.tracks = "generalTracks"
    
    process.slimmedJetsPFImpactParameterTagInfos = impactParameterTagInfos.clone()
    process.slimmedJetsPFImpactParameterTagInfos.jetTracks = "slimmedJetsPFJetTracksAssociatorAtVertex"
    process.slimmedJetsPFSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
    process.slimmedJetsPFSecondaryVertexTagInfos.trackIPTagInfos = "slimmedJetsPFImpactParameterTagInfos"
    
    process.slimmedJetsPFJetBtaggingSV = cms.Sequence(
    	process.slimmedJetsPFImpactParameterTagInfos *
    process.slimmedJetsPFSecondaryVertexTagInfos 
    )
    process.slimmedJetsPFJetsBtag = cms.Sequence(
    process.slimmedJetsPFJetTracksAssociatorAtVertex *
    process.slimmedJetsPFJetBtaggingSV
    )
    
    ## isotrack producer
    from SemiLeptonicWVA.Utils.trackIsolationMaker_cfi import trackIsolationFilter
    from SemiLeptonicWVA.Utils.trackIsolationMaker_cfi import trackIsolationCounter
    ## default
    process.IsolatedTracks = trackIsolationFilter.clone(
      doTrkIsoVeto= False,
      vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidatesTag = cms.InputTag("packedPFCandidates"),
      dR_ConeSize         = cms.double(0.3),
      dz_CutValue         = cms.double(0.05),
      minPt_PFCandidate   = cms.double(15.0),
      isoCut              = cms.double(0.1),
      )
    process.CountIsoTracks = trackIsolationCounter.clone(
      src = cms.InputTag("IsolatedTracks"),
      minNumber = 1,
      )

    process.substructureSequence = cms.Sequence()
    process.softdrop_onMiniAOD = cms.Sequence()
    process.pruning_onMiniAOD = cms.Sequence()
    process.redoPatJets = cms.Sequence()
    process.puppi_onMiniAOD = cms.Sequence()
    process.redoPuppiJets = cms.Sequence()

    if (doAK8Reclustering):
        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.chs = cms.EDFilter("CandPtrSelector",
                               src = cms.InputTag('packedPFCandidates'),
                               cut = cms.string('fromPV')
                               )

        process.NjettinessAK8 = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak8PFJetsCHS"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 0.8 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( -999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(-999),             # not used by default
                            akAxesR0 = cms.double(-999.0)        # not used by default
                            )

        process.ak4PFJetsCHS = ak4PFJetsCHS.clone(src = 'chs')
        process.ak8PFJetsCHS = ak8PFJetsCHS.clone( src = 'chs', jetPtMin = 100.0 )

        process.ak8PFJetsCHSPruned = ak8PFJetsCHSPruned.clone( src = 'chs', jetPtMin = 100.0 )
        process.ak8PFJetsCHSPrunedMass = ak8PFJetsCHSPrunedMass.clone()
        process.ak8PFJetsCHSSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'chs', jetPtMin = 100.0 )
        process.ak8PFJetsCHSSoftDropMass = ak8PFJetsCHSSoftDropMass.clone()

        process.substructureSequence+=process.chs
        process.substructureSequence+=process.ak8PFJetsCHS
        process.substructureSequence+=process.NjettinessAK8

        process.softdrop_onMiniAOD += process.ak8PFJetsCHSSoftDrop + process.ak8PFJetsCHSSoftDropMass
        process.pruning_onMiniAOD += process.ak8PFJetsCHSPruned + process.ak8PFJetsCHSPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from ak8PFJetsCHS

        process.patJetCorrFactorsAK8 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHS' )
        process.patJetsAK8 = patJetsAK8.clone( jetSource = 'ak8PFJetsCHS' )
        process.patJetsAK8.userData.userFloats.src = [ cms.InputTag("ak8PFJetsCHSPrunedMass"), cms.InputTag("ak8PFJetsCHSSoftDropMass"), cms.InputTag("NjettinessAK8:tau1"), cms.InputTag("NjettinessAK8:tau2"), cms.InputTag("NjettinessAK8:tau3")]
        process.patJetsAK8.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8") )
        process.selectedPatJetsAK8 = selectedPatJetsAK8.clone( cut = cms.string('pt > 20') )

        process.redoPatJets+=process.patJetCorrFactorsAK8
        process.redoPatJets+=process.patJetsAK8
        process.redoPatJets+=process.selectedPatJetsAK8

        if (reDoPruningAndSoftdrop):
            process.patJetCorrFactorsAK8Pruned = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHSPruned' )
            process.patJetsAK8Pruned = patJetsAK8.clone( jetSource = 'ak8PFJetsCHSPruned' )
            process.patJetsAK8Pruned.userData.userFloats.src = [ "" ]
            process.patJetsAK8Pruned.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8Pruned") )
            process.selectedPatJetsAK8Pruned = selectedPatJetsAK8.clone( 
                src = cms.InputTag('patJetsAK8Pruned'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsAK8Pruned
            process.redoPatJets+=process.patJetsAK8Pruned
            process.redoPatJets+=process.selectedPatJetsAK8Pruned

            process.patJetCorrFactorsAK8Softdrop = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHSSoftDrop' )
            process.patJetsAK8Softdrop = patJetsAK8.clone( jetSource = 'ak8PFJetsCHSSoftDrop' )
            process.patJetsAK8Softdrop.userData.userFloats.src = [ "" ]
            process.patJetsAK8Softdrop.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8Softdrop") )
            process.selectedPatJetsAK8Softdrop = selectedPatJetsAK8.clone(
                src = cms.InputTag('patJetsAK8Softdrop'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsAK8Softdrop
            process.redoPatJets+=process.patJetsAK8Softdrop
            process.redoPatJets+=process.selectedPatJetsAK8Softdrop

    if (doPuppi):

        from CommonTools.PileupAlgos.Puppi_cff import puppi
        from RecoJets.JetProducers.ak4PFJetsPuppi_cfi import ak4PFJetsPuppi

        process.ak8PFJetsPuppi = ak4PFJetsPuppi.clone( rParam = 0.8 )
        process.puppi = puppi.clone( candName = cms.InputTag('packedPFCandidates'),
                             vertexName = cms.InputTag('offlineSlimmedPrimaryVertices'))

        process.puppi_onMiniAOD = cms.Sequence(process.puppi + process.ak8PFJetsPuppi)

        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.NjettinessAK8 = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak8PFJetsPuppi"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 0.8 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( -999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(-999),             # not used by default
                            akAxesR0 = cms.double(-999.0)        # not used by default
                            )

        process.ak4PFJetsPuppi = ak4PFJetsPuppi.clone(src = 'puppi')
        process.ak8PFJetsPuppi = process.ak8PFJetsPuppi.clone( src = 'puppi', jetPtMin = 100.0 )

        process.ak8PFJetsCHSPruned = ak8PFJetsCHSPruned.clone( src = 'puppi', jetPtMin = 100.0 )
        process.ak8PFJetsCHSPrunedMass = ak8PFJetsCHSPrunedMass.clone(    
            matched = cms.InputTag("ak8PFJetsCHSPruned"),
            src = cms.InputTag("ak8PFJetsPuppi")
            )
        process.ak8PFJetsCHSSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'puppi', jetPtMin = 100.0 )
        process.ak8PFJetsCHSSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(
            matched = cms.InputTag("ak8PFJetsCHSSoftDrop"),
            src = cms.InputTag("ak8PFJetsPuppi")
            )

        process.substructureSequence+=process.puppi
        process.substructureSequence+=process.ak8PFJetsPuppi
        process.substructureSequence+=process.NjettinessAK8

        process.softdrop_onMiniAOD += process.ak8PFJetsCHSSoftDrop + process.ak8PFJetsCHSSoftDropMass
        process.pruning_onMiniAOD += process.ak8PFJetsCHSPruned + process.ak8PFJetsCHSPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from puppi AK8

        process.puppiJetCorrFactorsAK8 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsPuppi',
                                                                     levels = cms.vstring('L2Relative',
                                                                                          'L3Absolute')
                                                                     )
        process.puppiJetsAK8 = patJetsAK8.clone( jetSource = 'ak8PFJetsPuppi' )
        process.puppiJetsAK8.userData.userFloats.src = [ cms.InputTag("ak8PFJetsCHSPrunedMass"), cms.InputTag("ak8PFJetsCHSSoftDropMass"), cms.InputTag("NjettinessAK8:tau1"), cms.InputTag("NjettinessAK8:tau2"), cms.InputTag("NjettinessAK8:tau3")]
        process.puppiJetsAK8.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("puppiJetCorrFactorsAK8") )
        process.selectedPuppiJetsAK8 = selectedPatJetsAK8.clone( src = 'puppiJetsAK8', cut = cms.string('pt > 20') )

        process.redoPuppiJets+=process.puppiJetCorrFactorsAK8
        process.redoPuppiJets+=process.puppiJetsAK8
        process.redoPuppiJets+=process.selectedPuppiJetsAK8

#######AK8 GEN JETS################

    process.substructureSequenceGen = cms.Sequence()
    process.softdropGen_onMiniAOD = cms.Sequence()
    process.pruningGen_onMiniAOD = cms.Sequence()
    process.redoGenJets = cms.Sequence()
#    process.puppi_onMiniAOD = cms.Sequence()

    if (genJetsAK8Reclustering and MC):    
        from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
        
        process.ak8GenJets = ak4GenJets.clone(src = cms.InputTag('packedGenParticles'),
                                          rParam = cms.double(0.8)
                                          )

        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.NjettinessGenAK8 = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak8GenJets"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 0.8 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( -999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(-999),             # not used by default
                            akAxesR0 = cms.double(-999.0)        # not used by default
                            )

        process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
                                             src = cms.InputTag("packedGenParticles"),
                                             ignoreParticleIDs = cms.vuint32(
                1000022,
                1000012, 1000014, 1000016,
                2000012, 2000014, 2000016,
                1000039, 5100039,
                4000012, 4000014, 4000016,
                9900012, 9900014, 9900016,
                39),
                                             partonicFinalState = cms.bool(False),
                                             excludeResonances = cms.bool(False),
                                             excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
                                             tausAsJets = cms.bool(False)
                                             )

        from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters

        process.ak8GenJetsPruned = ak4GenJets.clone(
            SubJetParameters,
            rParam = cms.double(0.8),
            src = cms.InputTag("genParticlesForJets"),
            usePruning = cms.bool(True),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets")
            )

#        process.ak8GenJetsPruned = ak8PFJetsCHSPruned.clone( src = 'packedGenParticles', jetPtMin = 100.0 )
        process.ak8GenJetsPrunedMass = ak8PFJetsCHSPrunedMass.clone(    
            matched = cms.InputTag("ak8GenJetsPruned"),
            src = cms.InputTag("ak8GenJets")
            )

        process.ak8GenJetsSoftDrop = ak4GenJets.clone(
            SubJetParameters,
            rParam = cms.double(0.8),
            src = cms.InputTag("genParticlesForJets"),
            useSoftDrop = cms.bool(True),
            R0 = cms.double(0.8),
            beta = cms.double(0.0),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets")
            )

        process.ak8GenJetsSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(
            matched = cms.InputTag("ak8GenJetsSoftDrop"),
            src = cms.InputTag("ak8GenJets")
            )

        process.substructureSequenceGen+=process.genParticlesForJets
        process.substructureSequenceGen+=process.ak8GenJets
        process.substructureSequenceGen+=process.NjettinessGenAK8

        process.softdropGen_onMiniAOD += process.ak8GenJetsSoftDrop + process.ak8GenJetsSoftDropMass
        process.pruningGen_onMiniAOD += process.ak8GenJetsPruned + process.ak8GenJetsPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from gen AK8

        process.genJetsAK8 = patJetsAK8.clone( jetSource = 'ak8GenJets' )
        process.genJetsAK8.userData.userFloats.src = [ cms.InputTag("ak8GenJetsPrunedMass"), cms.InputTag("ak8GenJetsSoftDropMass"), cms.InputTag("NjettinessGenAK8:tau1"), cms.InputTag("NjettinessGenAK8:tau2"), cms.InputTag("NjettinessGenAK8:tau3")]
        process.genJetsAK8.addJetCorrFactors = cms.bool(False)
        process.genJetsAK8.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("") )
        process.selectedGenJetsAK8 = selectedPatJetsAK8.clone( src = 'genJetsAK8', cut = cms.string('pt > 20') )

        process.redoGenJets+=process.genJetsAK8
        process.redoGenJets+=process.selectedGenJetsAK8


######### A4PF-nonCHS jets ###########

    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

    process.ak4PFJets = ak4PFJets.clone(src = "packedPFCandidates")

    from SemiLeptonicWVA.Utils.ak4pfjets_cfi import patJetCorrFactorsAK4, patJetsAK4

    process.patJetCorrFactorsAK4 = patJetCorrFactorsAK4.clone( src = 'ak4PFJets' )
    process.patJetsAK4 = patJetsAK4.clone( jetSource = 'ak4PFJets' )


    #
    # Set up electron ID (VID framework)
    #
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
    # turn on VID producer, indicate data format to be
    # DataFormat.AOD or DataFormat.MiniAOD, as appropriate
    dataFormat=DataFormat.MiniAOD
    switchOnVIDElectronIdProducer(process,dataFormat)

    process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

    # define which IDs we want to produce
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
    #add them to the VID producer
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

    # Producers
    from SemiLeptonicWVA.Utils.electron_cfi import electron
    process.Electrons = electron.clone(
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        EleTag = cms.InputTag("slimmedElectrons"),
        MinPt = cms.double(-1),
        RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
        eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
        eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
        eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60")
    )

    # Add in Photons:
    switchOnVIDPhotonIdProducer(process, dataFormat)
    my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff']
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
    from SemiLeptonicWVA.Utils.photon_cfi import photon
    process.Photons = photon.clone(
        PhoTag = cms.InputTag("slimmedPhotons"),
        phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
        phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
        phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
        phoTightIdFullInfoMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight")
    )

    from SemiLeptonicWVA.Utils.muon_cfi import muon
    process.Muons = muon.clone(
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        MuTag = cms.InputTag("slimmedMuons"),
        MinPt = cms.double(-1),
        RhoTag = cms.InputTag("fixedGridRhoFastjetAll")
    )
    from SemiLeptonicWVA.Utils.subJetSelection_cfi import SubJetSelection
    process.HTJets = SubJetSelection.clone(
    JetTag  = cms.InputTag('slimmedJets'),
    MinPt								  = cms.double(50),
    MaxEta								  = cms.double(2.5),
    )
    from SemiLeptonicWVA.Utils.htdouble_cfi import htdouble
    process.HT = htdouble.clone(
    JetTag  = cms.InputTag('HTJets'),
    )
    from SemiLeptonicWVA.Utils.njetint_cfi import njetint
    process.NJets = njetint.clone(
    JetTag  = cms.InputTag('HTJets'),
    )
    from SemiLeptonicWVA.Utils.btagint_cfi import btagint
    process.BTags = btagint.clone(
    JetTag  = cms.InputTag('HTJets'),
    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    BTagCutValue					= cms.double(0.679)
    )
    from SemiLeptonicWVA.Utils.subJetSelection_cfi import SubJetSelection
    process.MHTJets = SubJetSelection.clone(
    JetTag  = cms.InputTag('slimmedJets'),
    MinPt								  = cms.double(30),
    MaxEta								  = cms.double(5.0),
    )
    process.MHTJetsAK8 = SubJetSelection.clone(
    JetTag  = cms.InputTag('slimmedJetsAK8'),
    MinPt								  = cms.double(30),
    MaxEta								  = cms.double(5.0),
    )
    from SemiLeptonicWVA.Utils.jetproperties_cfi import jetproperties
    process.MHTJetsProperties = jetproperties.clone(
    JetTag  = cms.InputTag('MHTJets'),
    doJEC  = cms.bool(doJECCorrection),
    L1File = cms.string("Summer15_25nsV7_DATA_L1FastJet_AK4PFchs.txt"),
    L2File = cms.string("Summer15_25nsV7_DATA_L2Relative_AK4PFchs.txt"),
    L3File = cms.string("Summer15_25nsV7_DATA_L3Absolute_AK4PFchs.txt"),
    L2L3File = cms.string("Summer15_25nsV7_DATA_L2L3Residual_AK4PFchs.txt"),
    )
    from SemiLeptonicWVA.Utils.jetpropertiesAK8_cfi import jetpropertiesAK8
    process.MHTJetsPropertiesAK8 = jetpropertiesAK8.clone(
    JetTag  = cms.InputTag('MHTJetsAK8'),
    puppiJetTag = cms.InputTag('selectedPuppiJetsAK8'),
    doJEC  = cms.bool(doJECCorrection),
    doReclusteringForPrunedAndSoftdrop = cms.bool(reDoPruningAndSoftdrop),
    L1File = cms.string("Summer15_25nsV7_DATA_L1FastJet_AK8PFchs.txt"),
    L2File = cms.string("Summer15_25nsV7_DATA_L2Relative_AK8PFchs.txt"),
    L3File = cms.string("Summer15_25nsV7_DATA_L3Absolute_AK8PFchs.txt"),
    L2L3File = cms.string("Summer15_25nsV7_DATA_L2L3Residual_AK8PFchs.txt"),
    )
    if (reDoPruningAndSoftdrop):
        process.MHTJetsPropertiesAK8.prunedJetTag  = cms.InputTag('selectedPatJetsAK8Pruned')
        process.MHTJetsPropertiesAK8.softdropJetTag  = cms.InputTag('selectedPatJetsAK8Softdrop')
    else:
        process.MHTJetsPropertiesAK8.prunedJetTag  = cms.InputTag('slimmedJetsAK8')
        process.MHTJetsPropertiesAK8.softdropJetTag  = cms.InputTag('slimmedJetsAK8')

    if (MC):
        process.MHTJetsProperties.L1File = cms.string("Summer15_25nsV7_MC_L1FastJet_AK4PFchs.txt")
        process.MHTJetsProperties.L2File = cms.string("Summer15_25nsV7_MC_L2Relative_AK4PFchs.txt")
        process.MHTJetsProperties.L3File = cms.string("Summer15_25nsV7_MC_L3Absolute_AK4PFchs.txt")
        process.MHTJetsProperties.L2L3File = cms.string("NONE")
        process.MHTJetsPropertiesAK8.L1File = cms.string("Summer15_25nsV7_MC_L1FastJet_AK8PFchs.txt")
        process.MHTJetsPropertiesAK8.L2File = cms.string("Summer15_25nsV7_MC_L2Relative_AK8PFchs.txt")
        process.MHTJetsPropertiesAK8.L3File = cms.string("Summer15_25nsV7_MC_L3Absolute_AK8PFchs.txt")
        process.MHTJetsPropertiesAK8.L2L3File = cms.string("NONE")

    from SemiLeptonicWVA.Utils.jetproperties_cfi import jetproperties
    process.JetsProperties = jetproperties.clone(
    JetTag  = cms.InputTag('slimmedJets'),
    MinPt = cms.double(-1),
    doJEC  = cms.bool(doJECCorrection),
    L1File = cms.string("Summer15_25nsV7_DATA_L1FastJet_AK4PFchs.txt"),
    L2File = cms.string("Summer15_25nsV7_DATA_L2Relative_AK4PFchs.txt"),
    L3File = cms.string("Summer15_25nsV7_DATA_L3Absolute_AK4PFchs.txt"),
    L2L3File = cms.string("Summer15_25nsV7_DATA_L2L3Residual_AK4PFchs.txt"),
    )
    from SemiLeptonicWVA.Utils.jetpropertiesAK8_cfi import jetpropertiesAK8
    process.JetsPropertiesAK8 = jetpropertiesAK8.clone(
    JetTag  = cms.InputTag('slimmedJetsAK8'),
    MinPt = cms.double(-1),
    doJEC  = cms.bool(doJECCorrection),
    doReclusteringForPrunedAndSoftdrop = cms.bool(reDoPruningAndSoftdrop),
    L1File = cms.string("Summer15_25nsV7_DATA_L1FastJet_AK8PFchs.txt"),
    L2File = cms.string("Summer15_25nsV7_DATA_L2Relative_AK8PFchs.txt"),
    L3File = cms.string("Summer15_25nsV7_DATA_L3Absolute_AK8PFchs.txt"),
    L2L3File = cms.string("Summer15_25nsV7_DATA_L2L3Residual_AK8PFchs.txt"),
    )
    if (reDoPruningAndSoftdrop):
        process.JetsPropertiesAK8.prunedJetTag  = cms.InputTag('selectedPatJetsAK8Pruned')
        process.JetsPropertiesAK8.softdropJetTag  = cms.InputTag('selectedPatJetsAK8Softdrop')
    else:
        process.JetsPropertiesAK8.prunedJetTag  = cms.InputTag('slimmedJetsAK8')
        process.JetsPropertiesAK8.softdropJetTag  = cms.InputTag('slimmedJetsAK8')
    if (MC):
        process.JetsProperties.L1File = cms.string("Summer15_25nsV7_MC_L1FastJet_AK4PFchs.txt")
        process.JetsProperties.L2File = cms.string("Summer15_25nsV7_MC_L2Relative_AK4PFchs.txt")
        process.JetsProperties.L3File = cms.string("Summer15_25nsV7_MC_L3Absolute_AK4PFchs.txt")
        process.JetsProperties.L2L3File = cms.string("NONE")
        process.JetsPropertiesAK8.L1File = cms.string("Summer15_25nsV7_MC_L1FastJet_AK8PFchs.txt")
        process.JetsPropertiesAK8.L2File = cms.string("Summer15_25nsV7_MC_L2Relative_AK8PFchs.txt")
        process.JetsPropertiesAK8.L3File = cms.string("Summer15_25nsV7_MC_L3Absolute_AK8PFchs.txt")
        process.JetsPropertiesAK8.L2L3File = cms.string("NONE")

    if doAK8Reclustering:
        process.JetsPropertiesAK8.JetTag = cms.InputTag('selectedPatJetsAK8')
    if doPuppi:
        process.JetsPropertiesAK8.JetTag = cms.InputTag('selectedPuppiJetsAK8')
    from SemiLeptonicWVA.Utils.mhtdouble_cfi import mhtdouble
    process.MHT = mhtdouble.clone(
    JetTag  = cms.InputTag('MHTJets'),
    )
    from SemiLeptonicWVA.Utils.deltaphidouble_cfi import deltaphidouble
    process.DeltaPhi = deltaphidouble.clone(
    DeltaPhiJets  = cms.InputTag('HTJets'),
    MHTJets  = cms.InputTag("MHTJets"),
    )
    from SemiLeptonicWVA.Utils.metdouble_cfi import metdouble
    process.MET = metdouble.clone(
    METTag  = cms.InputTag("slimmedMETs"),
    JetTag  = cms.InputTag('slimmedJets'),
    doJEC  = cms.bool(doJECCorrection),
    L1File = cms.string("Summer15_25nsV7_DATA_L1FastJet_AK4PFchs.txt"),
    L2File = cms.string("Summer15_25nsV7_DATA_L2Relative_AK4PFchs.txt"),
    L3File = cms.string("Summer15_25nsV7_DATA_L3Absolute_AK4PFchs.txt"),
    L2L3File = cms.string("Summer15_25nsV7_DATA_L2L3Residual_AK4PFchs.txt"),
    MuTag = cms.InputTag("slimmedMuons"),
    RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    corrMet = cms.bool(doJECCorrection),
    )

    if (MC):
        process.MET.L1File = cms.string("Summer15_25nsV7_MC_L1FastJet_AK4PFchs.txt")
        process.MET.L2File = cms.string("Summer15_25nsV7_MC_L2Relative_AK4PFchs.txt")
        process.MET.L3File = cms.string("Summer15_25nsV7_MC_L3Absolute_AK4PFchs.txt")
        process.MET.L2L3File = cms.string("NONE")

    from SemiLeptonicWVA.Utils.primaryverticies_cfi import primaryverticies
    process.NVtx = primaryverticies.clone(
    VertexCollection  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    )
    from SemiLeptonicWVA.Utils.genLeptonRecoCand_cfi import genLeptonRecoCand
    process.GenLeptons = genLeptonRecoCand.clone(
    PrunedGenParticleTag  = cms.InputTag("prunedGenParticles"),
    )
    from SemiLeptonicWVA.Utils.genJet_cfi import genJet
    process.GenJets = genJet.clone(
                            GenJetCollTag  = cms.InputTag("slimmedGenJets"),
                            )
    from SemiLeptonicWVA.Utils.genJetAK8_cfi import genJetAK8
    process.GenJetsAK8 = genJetAK8.clone(
                            GenJetCollTag  = cms.InputTag("selectedGenJetsAK8"),
                            )
    if not MC:
        process.GenLeptons = cms.Sequence()
        process.GenJets = cms.Sequence()
        process.GenJetsAK8 = cms.Sequence()



    ##### MET filters #####

    #### -----> MET Filter Flags from MiniAOD/TWiki <----- ####
    import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
    process.metBits_miniAOD = hlt.triggerResultsFilter.clone()
    # default is to use the latest process (but can set different process through Commandline Args)
    process.metBits_miniAOD.hltResults = cms.InputTag('TriggerResults::%s'%METFiltersProcess) 
    process.metBits_miniAOD.l1tResults = cms.InputTag('')
    #currently configured for CSCTightHaloFilter + GoodVertices
    met_bits = ['(Flag_CSCTightHaloFilter)','(Flag_goodVertices)','(Flag_eeBadScFilter)']
    bitsexpr = ' AND '.join(met_bits)
    process.metBits_miniAOD.triggerConditions = cms.vstring(bitsexpr)

    #### -----> HBHE noise filter <----- ####
    process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
    process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
    process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
    
    ########## save flags for filters
    from SemiLeptonicWVA.Utils.filterproducer_cfi import filterProducer
    process.FilterProducer = filterProducer.clone(
                                 noiseFilterTag = cms.InputTag("TriggerResults"),
                                 HBHENoiseFilter_Selector_ = cms.string("Flag_HBHENoiseFilter"),
                                 HBHENoiseIsoFilter_Selector_ = cms.string("Flag_HBHENoiseIsoFilter"),
                                 CSCHaloNoiseFilter_Selector_ = cms.string("Flag_CSCTightHaloFilter"),
                                 GoodVtxNoiseFilter_Selector_ = cms.string("Flag_goodVertices"),
                                 EEBadScNoiseFilter_Selector_ = cms.string("Flag_eeBadScFilter"),
                                 HBHENoiseFilterLoose = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Loose"),
                                 HBHENoiseFilterTight = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Tight"),
                                 HBHENoiseIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult")
                                 )


    RecoCandVector = cms.vstring()
    RecoCandVector.extend(['IsolatedTracks']) # basic muons electrons and isoalted tracks
    RecoCandVector.extend(['GenLeptons:Boson(GenBoson)|GenLeptons:BosonPDGId(I_GenBosonPDGId)','GenLeptons:Muon(GenMu)|GenLeptons:MuonTauDecay(I_GenMuFromTau)' ,'GenLeptons:Electron(GenElec)|GenLeptons:ElectronTauDecay(I_GenElecFromTau)','GenLeptons:Tau(GenTau)|GenLeptons:TauHadronic(I_GenTauHad)','GenLeptons:Neutrino(GenNu)'] ) # gen information on leptons
    RecoCandVector.extend(['GenJets:GenJet(GenJets)'] ) # gen information on jets
    RecoCandVector.extend(['GenJetsAK8:GenJetAK8(GenJetsAK8)|GenJetsAK8:GenAK8prunedMass(F_prunedMass)|GenJetsAK8:GenAK8softdropMass(F_softdropMass)|GenJetsAK8:GenAK8tau1(F_tau1)|GenJetsAK8:GenAK8tau2(F_tau2)|GenJetsAK8:GenAK8tau3(F_tau3)'] ) # gen information on AK8 jets
    RecoCandVector.extend(['JetsProperties(Jets)|JetsProperties:bDiscriminatorCSV(F_bDiscriminatorCSV)|JetsProperties:bDiscriminatorICSV(F_bDiscriminatorICSV)|JetsProperties:chargedEmEnergyFraction(F_chargedEmEnergyFraction)|JetsProperties:chargedHadronEnergyFraction(F_chargedHadronEnergyFraction)|JetsProperties:chargedHadronMultiplicity(I_chargedHadronMultiplicity)|JetsProperties:electronMultiplicity(I_electronMultiplicity)|JetsProperties:jetArea(F_jetArea)|JetsProperties:muonEnergyFraction(F_muonEnergyFraction)|JetsProperties:muonMultiplicity(I_muonMultiplicity)|JetsProperties:neutralEmEnergyFraction(F_neutralEmEnergyFraction)|JetsProperties:neutralHadronMultiplicity(I_neutralHadronMultiplicity)|JetsProperties:photonEnergyFraction(F_photonEnergyFraction)|JetsProperties:photonMultiplicity(I_photonMultiplicity)|JetsProperties:isLooseJetId(b_isLooseJetId)|JetsProperties:isTightJetId(b_isTightJetId)|JetsProperties:isTightLepVetoJetId(b_isTightLepVetoJetId)|JetsProperties:PtCorr(F_PtCorr)|JetsProperties:EtaCorr(F_EtaCorr)|JetsProperties:PhiCorr(F_PhiCorr)|JetsProperties:ECorr(F_ECorr)'] ) # jet information on various variables
    RecoCandVector.extend(['JetsPropertiesAK8(AK8Jets)|JetsPropertiesAK8:AK8bDiscriminatorCSV(F_bDiscriminatorCSV)|JetsPropertiesAK8:AK8bDiscriminatorICSV(F_bDiscriminatorICSV)|JetsPropertiesAK8:AK8chargedEmEnergyFraction(F_chargedEmEnergyFraction)|JetsPropertiesAK8:AK8chargedHadronEnergyFraction(F_chargedHadronEnergyFraction)|JetsPropertiesAK8:AK8chargedHadronMultiplicity(I_chargedHadronMultiplicity)|JetsPropertiesAK8:AK8electronMultiplicity(I_electronMultiplicity)|JetsPropertiesAK8:AK8jetArea(F_jetArea)|JetsPropertiesAK8:AK8muonEnergyFraction(F_muonEnergyFraction)|JetsPropertiesAK8:AK8muonMultiplicity(I_muonMultiplicity)|JetsPropertiesAK8:AK8neutralEmEnergyFraction(F_neutralEmEnergyFraction)|JetsPropertiesAK8:AK8neutralHadronMultiplicity(I_neutralHadronMultiplicity)|JetsPropertiesAK8:AK8photonEnergyFraction(F_photonEnergyFraction)|JetsPropertiesAK8:AK8photonMultiplicity(I_photonMultiplicity)|JetsPropertiesAK8:AK8prunedMass(F_prunedMass)|JetsPropertiesAK8:AK8softDropMass(F_softDropMass)|JetsPropertiesAK8:AK8trimmedMass(F_trimmedMass)|JetsPropertiesAK8:AK8filteredMass(F_filteredMass)|JetsPropertiesAK8:AK8tau1(F_tau1)|JetsPropertiesAK8:AK8tau2(F_tau2)|JetsPropertiesAK8:AK8tau3(F_tau3)|JetsPropertiesAK8:AK8isLooseJetId(b_AK8isLooseJetId)|JetsPropertiesAK8:AK8isTightJetId(b_AK8isTightJetId)|JetsPropertiesAK8:AK8isTightLepVetoJetId(b_AK8isTightLepVetoJetId)|JetsPropertiesAK8:PtCorr(F_PtCorr)|JetsPropertiesAK8:EtaCorr(F_EtaCorr)|JetsPropertiesAK8:PhiCorr(F_PhiCorr)|JetsPropertiesAK8:ECorr(F_ECorr)'] ) # AK8 jet information on various variables
# 
    RecoCandVector.extend(['Electrons(Electrons)|Electrons:charge(I_charge)|Electrons:isHEEP(b_isHEEP)|Electrons:type(I_type)|Electrons:mass(F_mass)|Electrons:pfDeltaCorrRelIso(F_pfDeltaCorrRelIso)|Electrons:pfRhoCorrRelIso04(F_pfRhoCorrRelIso04)|Electrons:pfRhoCorrRelIso03(F_pfRhoCorrRelIso03)|Electrons:pfRelIso(F_pfRelIso)|Electrons:photonIso(F_photonIso)|Electrons:neutralHadIso(F_neutralHadIso)|Electrons:chargedHadIso(F_chargedHadIso)|Electrons:trackIso(F_trackIso)|Electrons:isLoose(b_isLoose)|Electrons:isMedium(b_isMedium)|Electrons:isTight(b_isTight)|Electrons:SCEnergy(F_SCEnergy)|Electrons:deltaEtaSCTracker(F_deltaEtaSCTracker)|Electrons:deltaPhiSCTracker(F_deltaPhiSCTracker)|Electrons:sigmaIetaIeta(F_sigmaIetaIeta)|Electrons:sigmaIphiIphi(F_sigmaIphiIphi)'] ) # electron information on various variables
    RecoCandVector.extend(['Muons(Muons)|Muons:charge(I_charge)|Muons:isHighPt(b_isHighPt)|Muons:type(I_type)|Muons:mass(F_mass)|Muons:pfDeltaCorrRelIso(F_pfDeltaCorrRelIso)|Muons:pfRelIso(F_pfRelIso)|Muons:photonIso(F_photonIso)|Muons:neutralHadIso(F_neutralHadIso)|Muons:chargedHadIso(F_chargedHadIso)|Muons:trackIso(F_trackIso)|Muons:isLoose(b_isLoose)|Muons:isMedium(b_isMedium)|Muons:isTight(b_isTight)|Muons:isPFMuon(b_isPFMuon)'] ) # muon information on various variables
    RecoCandVector.extend(['Photons(Photons)|Photons:isLoose(b_isLoose)|Photons:isMedium(b_isMedium)|Photons:isTight(b_isTight)|Photons:minPt(F_minPt)|Photons:phoSCEtaMultiRange(F_phoSCEtaMultiRange)|Photons:phoSingleTowerHadOverEm(F_phoSingleTowerHadOverEm)|Photons:phoFull5x5SigmaIEtaIEta(F_phoFull5x5SigmaIEtaIEta)|Photons:phoAnyPFIsoWithEA(F_phoAnyPFIsoWithEA)|Photons:phoAnyPFIsoWithEAAndExpoScaling(F_phoAnyPFIsoWithEAAndExpoScaling)|Photons:phoAnyPFIsoWithEA1(F_phoAnyPFIsoWithEA1)|Photons:hasPixelSeed(b_hasPixelSeed)|Photons:passElectronVeto(b_passElectronVeto)|Photons:photonIso(F_photonIso)|Photons:neutralHadIso(F_neutralHadIso)|Photons:chargedHadIso(F_chargedHadIso)|Photons:puChargedHadIso(F_puChargedHadIso)|Photons:sigmaIetaIeta(F_sigmaIetaIeta)'] ) # photon information on various variables

    from SemiLeptonicWVA.TreeMaker.treeMaker import TreeMaker
    process.TreeMaker2 = TreeMaker.clone(
    	TreeName          = cms.string("PreSelection"),
    	VarsRecoCand = RecoCandVector,
    	VarsDouble  	  = cms.vstring('WeightProducer:weight(Weight)','MHT','MET:Pt(METPt)','MET:Phi(METPhi)','MET:PtRaw(METPtRaw)','MET:PhiRaw(METPhiRaw)','MET:CaloMetPt(CaloMetPt)','MET:CaloMetPhi(CaloMetPhi)','HT','DeltaPhi:DeltaPhi1(DeltaPhi1)','DeltaPhi:DeltaPhi2(DeltaPhi2)','DeltaPhi:DeltaPhi3(DeltaPhi3)','GenEventInfo:genEventWeight(genEventWeight)','GenEventInfo:PUWeight(PUWeight)','GenEventInfo:originalWeight(originalWeight)'),
    	VarsInt = cms.vstring('NJets','BTags','NVtx','GenEventInfo:npT(npT)','FilterProducer:passFilterHBHE(passFilterHBHE)','FilterProducer:passFilterHBHEIso(passFilterHBHEIso)','FilterProducer:passFilterCSCHalo(passFilterCSCHalo)','FilterProducer:passFilterGoodVtx(passFilterGoodVtx)','FilterProducer:passFilterEEBadSC(passFilterEEBadSC)','FilterProducer:passFilterHBHELooseRerun(passFilterHBHELooseRerun)','FilterProducer:passFilterHBHETightRerun(passFilterHBHETightRerun)','FilterProducer:passFilterHBHEIsoRerun(passFilterHBHEIsoRerun)'),
        debug = debug,
    	)

    process.TreeMaker2.VarsDouble.extend(['MET:PtDefault(METPtDefault)','MET:PhiDefault(METPhiDefault)','MET:PtType1(METPtType1)','MET:PhiType1(METPhiType1)','MET:PtType1XYSmear(METPtType1XYSmear)','MET:PhiType1XYSmear(METPhiType1XYSmear)','MET:PtType1Smear(METPtType1Smear)','MET:PhiType1Smear(METPhiType1Smear)','MET:PtType1XY(METPtType1XY)','MET:PhiType1XY(METPhiType1XY)'])
    process.TreeMaker2.VectorBool.extend(['TriggerProducer:TriggerPass'])
    process.TreeMaker2.VectorInt.extend(['TriggerProducer:TriggerPrescales'])
    process.TreeMaker2.VectorDouble.extend(['GenEventInfo:AQGCweights'])
    process.TreeMaker2.VectorString.extend(['TriggerProducer:TriggerNames'])

    ## --- Final paths ----------------------------------------------------
    process.out = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string("output.root"),
                                   )
        
    process.dump = cms.EDAnalyzer("EventContentAnalyzer")
    process.WriteTree = cms.Path(
      process.TriggerProducer*
      ### MET Filter Bits
      process.HBHENoiseFilterResultProducer*
      process.FilterProducer* #this now contains all the met filters
      ### rest of ntupling starts after here
      process.filterSeq *
      process.GenEventInfo *
      process.Muons *
      process.egmGsfElectronIDSequence*
      process.Electrons *
      process.egmPhotonIDSequence *
      process.Photons *
      process.WeightProducer *
      process.IsolatedTracks *
      process.substructureSequenceGen *
      process.softdropGen_onMiniAOD *
      process.pruningGen_onMiniAOD *
      process.redoGenJets*
      process.GenJetsAK8 *
      process.puppi_onMiniAOD *
      process.substructureSequence *
      process.softdrop_onMiniAOD *
      process.pruning_onMiniAOD *
      process.redoPatJets*
      process.redoPuppiJets*
      process.HTJets *
      process.HT *
      process.NJets *
      process.BTags *
      process.MHTJets *
      process.MHTJetsProperties *
      process.JetsProperties *
      process.MHTJetsAK8 *
      process.MHTJetsPropertiesAK8 *
      process.JetsPropertiesAK8 *
      process.MHT *
      process.MET *
      process.DeltaPhi *
      process.NVtx *
      process.GenLeptons *
      process.GenJets *
      process.TreeMaker2

        )
