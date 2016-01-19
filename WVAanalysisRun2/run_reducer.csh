#!/bin/csh -f

set input_dir = "/eos/uscms/store/user/jfaulkn3/PostFailure_Files/"
set MC = ( TTGJets TGJets WGjets WGjets_PtG500 ZGjets WWG WZG )
set XSec = ( 3.697 2.967 405.271 0.0117887 117.864 0.2147 0.04123 )
set Nevts = ( 4843746 400000 6099599 1392843 4455517 1000000 1000000 )

foreach lep ( el )

  foreach mc ( `seq 1 $#MC` )

    set com = "python/produceWWNtuples.py -i $input_dir -n MC_25ns_3p8T_$MC[$mc].root -o RD_$MC[$mc]_${lep}.root -mc True -l $lep -no $Nevts[$mc] -w $XSec[$mc]"
    echo $com
    #python $com

  end

  if ($lep == "el") then
    set data_tag = Electron
  else
    set data_tag = Muon
  endif

  set com = "python/produceWWNtuples.py -i $input_dir -n Run2015CD_25ns_3p8T_Single${data_tag}.root -o RD_Run2015CD_Single${data_tag}.root -mc False -l $lep -no 10000000000"
  echo $com
  python $com

end

#EOF
