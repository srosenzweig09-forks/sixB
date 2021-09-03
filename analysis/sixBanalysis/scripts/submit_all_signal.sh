ODIR="/store/user/srosenzw/signal_exploration"

. scripts/arg_submit.sh -v sr "$@"
TAG="${TAG}NMSSM"

make exe -j || exit -1

echo "... tag       : ", $TAG
echo "... saving to : ", $ODIR

files=$(ls input/PrivateMC_2021/NMSSM_XYH_YToHH_6b_MX_700_MY_400.txt )
#input/PrivateMC_2021/NMSSM_XYH_YToHH_6b_MX_700_MY_400_10M.txt

for input in ${files[@]}; do
    python scripts/submitSkimOnBatch.py --tag $TAG --outputDir $ODIR --cfg $CFG --njobs 100 --input $input --is-signal
done
