
set -e
make exe -j || exit -1

relDir='scores/'
files=$(ls ${relDir})
for input in ${files[@]}; do
    infile="${input%.*}"

    if [[ "$infile" == "data" ]]; then
        # bin/gnn_jet_pairing.exe --is-data
        continue
    fi

    bin/gnn_jet_pairing.exe --is-signal --score-file $infile
done