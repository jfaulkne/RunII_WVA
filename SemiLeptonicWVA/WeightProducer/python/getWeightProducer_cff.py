# $Id: getWeightProducer_cff.py,v 1.7 2012/11/05 14:58:14 mschrode Exp $
#
# Returns a WeightProducer module that knows at runtime
# which data sample is produced and thus, what weights
# are required. The function can be used as follows:
#
#    from RA2Classic.WeightProducer.getWeightProducer_cff import getWeightProducer
#    process.WeightProducer = getWeightProducer(process.source.fileNames[0])


import FWCore.ParameterSet.Config as cms

def getWeightProducer(fileName):

    mcVersion = "none"                  # For lumi and PU weights
    applyWeight = False


    ## --- Setup default WeightProducer ------------------------------------

    # Import weightProducer
    from SemiLeptonicWVA.WeightProducer.weightProducer_cfi import weightProducer

    # Set default values to produce an event weight of 1
    weightProducer.weight = cms.double(1.0)
    weightProducer.Method = cms.string("")
    weightProducer.LumiScale = cms.double(1.0)
    weightProducer.FileNamePUDataDistribution = cms.string("NONE")
    weightProducer.PU = cms.int32(0)



    ## --- Adjust WeightProducer for MC samples --------------------
    if "WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" in fileName:
                        mcVersion = "Spring15MiniAODv2"
                        weightProducer.Method     = cms.string("Constant")
                        weightProducer.XS         = cms.double(405.271)
                        weightProducer.NumberEvts = cms.double(6099599)
                        applyWeight = True
                        weightProducer.weight = cms.double(-1.)

    if "WGToLNuG_PtG-500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" in fileName:
                        mcVersion = "Spring15MiniAODv2"
                        weightProducer.Method     = cms.string("Constant")
                        weightProducer.XS         = cms.double(0.0117887)
                        weightProducer.NumberEvts = cms.double(1391312)
                        applyWeight = True
                        weightProducer.weight = cms.double(-1.)

    if "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" in fileName:
                        mcVersion = "Spring15MiniAODv2"
                        weightProducer.Method     = cms.string("Constant")
                        weightProducer.XS         = cms.double(117.864)
                        weightProducer.NumberEvts = cms.double(4451319)
                        applyWeight = True
                        weightProducer.weight = cms.double(-1.)

    if "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8" in fileName:
                        mcVersion = "Spring15MiniAODv2"
                        weightProducer.Method     = cms.string("Constant")
                        weightProducer.XS         = cms.double(3.697)
                        weightProducer.NumberEvts = cms.double(4832230)
                        applyWeight = True
                        weightProducer.weight = cms.double(-1.)

    if "WZAToLNu2jA_4f_NLO_TuneCUETP8M1_13TeV-amcatnlo-pythia8" in fileName:
                        mcVersion = "Spring15MiniAODv2"
                        weightProducer.Method     = cms.string("Constant")
                        weightProducer.XS         = cms.double(0.04123)
                        weightProducer.NumberEvts = cms.double(997000)
                        applyWeight = True
                        weightProducer.weight = cms.double(-1.)


        
    ## --- PU Reweighting and Lumi ------------------------------------------------
         
    if mcVersion == "Summer12_5_3_X":
        weightProducer.weight = cms.double(-1.)
        weightProducer.Lumi = cms.double(19466)
        weightProducer.PU = cms.int32(3) # PU S10
## weightProducer.FileNamePUDataDistribution = cms.string("RA2Classic/AdditionalInputFiles/DataPileupHistogram_RA2Summer12_5_2_X_190456-196531_8TeV_PromptReco_WOLowPU.root")
        # grid-control requires additional input files to be in the data directory
        weightProducer.FileNamePUDataDistribution = cms.string("RA2Classic/WeightProducer/data/DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root")
        applyWeight = True
	
	## use this for PU SU7
#    if mcVersion == "Summer12_5_3_X":
#        weightProducer.weight = cms.double(-1.)
#        weightProducer.Lumi = cms.double(19466)
#        weightProducer.PU = cms.int32(2) # PU S7
## weightProducer.FileNamePUDataDistribution = cms.string("RA2Classic/AdditionalInputFiles/DataPileupHistogram_RA2Summer12_5_2_X_190456-196531_8TeV_PromptReco_WOLowPU.root")
        # grid-control requires additional input files to be in the data directory
        weightProducer.FileNamePUDataDistribution = cms.string("RA2Classic/WeightProducer/data/DataPUProfile_2013Jan22_XS69400_PixelLumiCorr_MediumMass.root")
        applyWeight = True
    if applyWeight:
        print "Setup WeightProducer for '"+fileName+"'"

    return weightProducer

    
