#!/bin/sh

# eval `scramv1 runtime -sh`
# source scripts/setup.sh

# output="NMSSM_XYH_YToHH_6b_MX_700_MY_400_accstudies_500k_Jul2021-v2.root"
# output="output-tree_unblind.root"
output="NMSSM_XYY_YToHH_8b_MX_1000_MY_450_accstudies.root"


# input="input/Run2_UL/2018/JetHT_Run2018C.txt --is-data "
# input="input/Run2_UL/2018/QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8.txt"
# input="input/PrivateMC_2021/NMSSM_XYH_YToHH_6b_MX_700_MY_400.txt --is-signal --no-genw-tree"
input="input/PrivateMC_2018/NMSSM_XYY_YToHH_8b/NMSSM_XYY_YToHH_8b_MX_1000_MY_450.txt --is-signal"
# input="input/PrivateMC_2018/NMSSM_XYY_YToHH_8b/training_5M/NMSSM_XYY_YToHH_8b_MX_1000_MY_450.txt --is-signal"

cfg="config/8b_config/skim_ntuple_2018_ranked_quadh.cfg"
# cfg="config/8b_config/skim_ntuple_2018_t8btag.cfg"
# cfg="config/8b_config/testing.cfg"
# cfg="config/8b_config/skim_ntuple_2018_gnn_jet.cfg"

# exe="perf record -g bin/skim_ntuple.exe"
exe=bin/skim_ntuple.exe

make exe -j && \
    $exe \
	--input $input \
	--cfg  $cfg \
	--output $output \
	$@
