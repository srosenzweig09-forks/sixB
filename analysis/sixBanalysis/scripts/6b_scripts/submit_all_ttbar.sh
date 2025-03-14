ODIR="/store/user/srosenzw/sixb/ntuples/Summer2018UL/"

# . scripts/arg_submit.sh -v qcd "$@"
TAG="dnn"
CFG="config/skim_ntuple_2018_marina.cfg"

make exe -j || exit -1

echo "... tag       : ", $TAG
echo "... saving to : ", $ODIR

python scripts/submitSkimOnBatch.py --tag $TAG --outputDir $ODIR --cfg $CFG --njobs 100 --input input/Run2_UL/2018/TTJets.txt --forceOverwrite   
# # python scripts/submitSkimOnBatch.py --tag $TAG --outputDir $ODIR --cfg config/skim_ntuple_2018_ttbar.cfg --njobs 100 --input input/Run2_UL/2018/SingleMuon_Run2.txt --is-data


#############

# TAG="ttbar_2018_10Jan2022"
# ODIR="/store/group/lpchbb/lcadamur/sixb_ntuples/"

# echo "... tag       : ", $TAG
# echo "... saving to : ", $ODIR

# python scripts/submitSkimOnBatch.py --tag $TAG --outputDir $ODIR --cfg config/skim_ntuple_2018_ttbar.cfg --njobs 100 --input input/Run2_UL/2018/TTJets.txt         
# python scripts/submitSkimOnBatch.py --tag $TAG --outputDir $ODIR --cfg config/skim_ntuple_2018_ttbar.cfg --njobs 100 --input input/Run2_UL/2018/SingleMuon_Run2.txt --is-data