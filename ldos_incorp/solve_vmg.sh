#!/bin/bash

cwd=$(pwd)
base_name="config"
num_cores=4
RATE_FILE=$(cat ${base_name}.yml | shyaml get-value rate_file)

if [ -z "$1" ]; then
    echo "No argument supplied, using run_dir=Sol/"
    run_dir=Sol # run_dir in this case is the twist angle folder with rate data
    if [ ! -d $run_dir ]; then
    echo "Run directory ${run_dir} does not exist. Making new dir Sol/"
    mkdir -p $run_dir
    fi
else
    run_dir=$1
    if [ ! -d $run_dir ]; then
    echo "Run directory ${run_dir} does not exist. Making new dir"
    mkdir -p $run_dir
    fi
fi

#$(cat ${base_name}.yml | shyaml get-value sol_dir)

# # Define Vapp values
# low=0.0
# high=0.0
# incr=0.1
# Vapp=( $(seq $low $incr $high) )


cd $run_dir
out=iv_data.txt
if [ -f "$out" ]; then
    echo "$out exists. Overwriting"
    echo "I Vapp" > $out
else
    echo "I Vapp" > $out
fi

for f in *.mat  #"${!Vapp[@]}"
do 
    fname=${f%.*}
    export SOL_DIR="${cwd}/${run_dir}"
    export SOL_FILE="${SOL_DIR}/config_${fname}"
    export RATE_LOC="${SOL_DIR}/${f}"
    FILE=${SOL_FILE}.h5
    if [ -f "$FILE" ]; then
        echo -e "${SOL_FILE} exists. Going to next config \n"
        continue
    else
        echo -e "Generating config file with rate data in ${SOL_DIR}/${f}"
        cp ${base_name}.yml ${SOL_FILE}.yml
        envsubst < ${cwd}/${base_name}.yml > ${SOL_FILE}.yml
        echo $SOL_FILE
        echo -e "Running ${fname} \n"
        #/Applications/Julia-1.6.app/Contents/Resources/julia/bin/julia $RATE_FILE $SOL_FILE
        mpirun -np $num_cores python ../pnp_par.py $SOL_FILE
        echo -e "Ended ${fname} \n"
    fi
done
cd $cwd