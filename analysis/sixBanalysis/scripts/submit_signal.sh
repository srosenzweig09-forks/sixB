
ODIR="/store/user/srosenzw/sixb_ntuples"
TAG="preselections"
input="fileList.txt"

make exe -j || exit -1

echo "... tag       : ", $TAG
echo "... saving to : ", $ODIR

python scripts/submitSkimOnBatch.py --tag $TAG --outputDir $ODIR --cfg config/skim_ntuple_2018.cfg --njobs 100 --input $input --is-signal --outputName NMSSM_XYH_YToHH_6b_MX_700_MY_400

# files=$(ls input/PrivateMC_2021/*)

# for input in ${files[@]}; do
    # python scripts/submitSkimOnBatch.py --tag $TAG --outputDir $ODIR --cfg config/skim_ntuple_2018.cfg --njobs 100 --input $input --is-signal
# done
