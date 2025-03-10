#!/bin/sh

# eval `scramv1 runtime -sh`
# source scripts/setup.sh

# output="NMSSM_XYH_YToHH_6b_MX_700_MY_400_accstudies_500k_Jul2021-v2.root"
# output="output-tree.root"
output="output.root"

# input="input/Run2_UL/2018/JetHT_Run2018C.txt --is-data "
# input="input/Run2_UL/2018/QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8.txt"
input="input/PrivateMC_2018/NMSSM_XYH_YToHH_6b/NMSSM_XYH_YToHH_6b_MX_700_MY_400.txt --is-signal"

cfg="config/skim_ntuple_2018_presel.cfg"
# cfg="config/skim_ntuple_2018_8b.cfg"
# cfg="config/skim_ntuple_2018_qcd.cfg"
# cfg="config/skim_ntuple_2018_cr.cfg"
     
make exe -j && \
    $exe \
	--input $input \
	--cfg  $cfg \
	--output $output \
	--maxEvts 100 \
	$@
