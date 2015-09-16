#!/bin/bash
export PYTHONHOME=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/python/2.7.6-cms

CHS=("el" "mu")
ERAS=("july" "prompt" "254833")
STAGES_EL=(24 43 13)
STAGES_MU=(24 27 9)

for CH in "${CHS[@]}"; do

   tick=1
   for ((index = 0; index < 3; index++)); do

      if [ "$CH" == "mu" ]; then
         STAGE=${STAGES_MU[$index]}
      else STAGE=${STAGES_EL[$index]}
      fi

      for ((run = 0; run < $STAGE; run++)); do

         bsub -q 8nh cmsRun /afs/cern.ch/work/j/jfaulkne/CMSSW_7_4_7_patch2/src/AllHadronicSUSY/TreeMaker/test/runMakeTreeFromMiniAOD_cfg.py channel=${CH} era=${ERAS[$index]} run=$run outfile=/afs/cern.ch/work/j/jfaulkne/CMSSW_7_4_7_patch2/src/AllHadronicSUSY/ReducedSelection_${tick}_${CH}
	 #cmsRun TreeMaker/test/runMakeTreeFromMiniAOD_cfg.py channel=${CH} era=${ERAS[$index]} run=$run outfile=ReducedSelection_${CH}_${tick}
         tick=$((tick+1))

      done

   done

done

#end
