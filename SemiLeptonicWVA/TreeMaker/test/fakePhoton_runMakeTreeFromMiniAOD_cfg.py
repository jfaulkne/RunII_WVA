# $Id: runMakeTreeFromPAT_cfg.py,v 1.9 2013/01/24 15:42:54 mschrode Exp $
#
# Expects a file name as argument e.g.
# cmsRun runMakeLostLeptonTreeFromPAT_cfg.py dataset=/store/user/mschrode/HT/RA2PreSelection_Run2012A-13Jul2012-v1_V4/21a074f94cdbe7cfbeeb19be46b40a6a/RA2Skim_9_1_h6A.root
# cmsRun ../test/runMakeLostLeptonTreeFromPAT_cfg.py dataset=/store/user/mschrode/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/RA2PreSelection_Summer12_DR53X-PU_S10_START53_V7A-v1_V4/6c50609e978ba7d5388d5439fc628605/RA2Skim_100_1_dgv.root, global_tag=START53_V7F::All, MC=True, debug=True

# Read parameters
from SemiLeptonicWVA.Utils.CommandLineParams import CommandLineParams
parameters = CommandLineParams()

MC= parameters.value("MC", False)
CH = parameters.value("channel","el")
RUN = parameters.value("run",0)
ERA = parameters.value("era","july")
OUTFILE = parameters.value("outfile","ReducedSelection")

if MC:
   #from SemiLeptonicWVA.Utils.MC_list import MC_LIST
   dataSetName = parameters.value("dataset", "")#MC_LIST[ERA][RUN])
else:
   #from SemiLeptonicWVA.Utils.DATA_list import SingleLepton_Dataset
   dataSetName = parameters.value("dataset", "")#SingleLepton_Dataset[CH][ERA][RUN])

JSON = {"prompt":"json/Cert_246908-254879_13TeV_PromptReco_Collisions15_JSON_prompt.txt",
"july":"json/Cert_246908-254879_13TeV_PromptReco_Collisions15_JSON_july.txt",
"254833":"json/Cert_254833_13TeV_PromptReco_Collisions15_JSON.txt"}

global_tag = parameters.value("global_tag","74X_dataRun2_v2")
QCD= parameters.value("QCD", False)
LostLepton= parameters.value("LostLepton", False)
debug= parameters.value("debug", False)
nJetsMin    = parameters.value("njets_min",0)
htMin       = parameters.value("ht_min",0)
mhtMin      = parameters.value("mht_min",0)
NumProcessedEvt=parameters.value("NumProcessedEvt",-1)
METFiltersProcess=parameters.value("METFiltersProcess","")
DoAK8Reclustering=parameters.value("DoAK8Reclustering",True)
DoJECCorrection=parameters.value("DoJECCorrection",True)
DoPuppi=parameters.value("DoPuppi",False)
LeptonFilter=parameters.value("leptonFilter",True)
GenJetsAK8Reclustering=parameters.value("genJetsAK8Reclustering",True)
isHBHEEarlyData = parameters.value("isHBHEEarlyData",True)
if MC: json = ""
else: json = JSON[ERA]
JsonFileName=parameters.value("jsonFileName",json)
IsCrab=parameters.value("isCrab",False)

processName      = parameters.value("name","WVGamma")


print "***** SETUP ************************************"
print "  dataSetName : "+dataSetName
print " global_tag : "+global_tag
print " runningOnMC : "+str(MC)
print " runningOnQCD : "+str(QCD)
print " LostLepton(MC) : "+str(LostLepton)
print "     nJetsMin : "+str(nJetsMin)
print "        htMin : "+str(htMin)
print "       mhtMin : "+str(mhtMin)
print "       debug : "+str(debug)
print "       num of events : "+str(NumProcessedEvt)
print "       doAK8Reclustering : "+str(DoAK8Reclustering)
print "       doJECCorrection : "+str(DoJECCorrection)
print "       doPuppi : "+str(DoPuppi)
print "       leptonFilter : "+str(LeptonFilter)
print "       genJetsAK8Reclustering : "+str(GenJetsAK8Reclustering)
print "       isHBHEEarlyData : "+str(isHBHEEarlyData)
print "       jsonFileName : "+str(JsonFileName)
print "       isCrab : "+str(False)
print "************************************************"

# The process needs to be defined AFTER reading sys.argv,
# otherwise edmConfigHash fails
import FWCore.ParameterSet.Config as cms
#process = cms.Process("RA2EventSelection")
process = cms.Process(processName)

from SemiLeptonicWVA.TreeMaker.fakePhoton_makeTreeFromMiniAOD_cff import makeTreeTreeFromMiniAOD
makeTreeTreeFromMiniAOD(process,
                outFileName=OUTFILE,
                NJetsMin=nJetsMin,
                HTMin=htMin,
                MHTMin=mhtMin,
                reportEveryEvt=1000,
                testFileName=dataSetName,
		Global_Tag=global_tag,
                METFiltersProcess=METFiltersProcess,
		MC=MC,
		QCD=QCD,
		LostLepton=LostLepton,
		debug = debug,
                numProcessedEvt=NumProcessedEvt,
                doAK8Reclustering=DoAK8Reclustering,
                doJECCorrection=DoJECCorrection,
                doPuppi=DoPuppi,
                leptonFilter=LeptonFilter,
                genJetsAK8Reclustering=GenJetsAK8Reclustering,
                customizeHBHENoiseForEarlyData=isHBHEEarlyData,
                jsonFileName=JsonFileName,
                isCrab=IsCrab)

