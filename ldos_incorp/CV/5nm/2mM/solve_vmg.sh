#!/bin/bash

start=`date +%s`

cwd=$(pwd)
base_name="config"
out=iv_data.txt
num_cores=6
RATE_FILE=$(cat ${base_name}.yml | shyaml get-value rate_file)

if [ -z "$1" ]; then
    echo "No argument supplied, using run_dir=Sol/"
    run_dir=Sol
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

##$(cat ${base_name}.yml | shyaml get-value sol_dir)
## Create out file wit IV data
cd $run_dir
FILE=iv_data.txt
if [ -f "$FILE" ]; then
    echo "$FILE exists. Overwriting"
    echo "Domain I Vapp" > $out
else
    echo "I Vapp Vdl" > $out
fi

## Defined Vapp values in run_jobs.py
# low=-0.4
# high=0.3
# incr=0.05
# Vapp=( $(seq $low $incr $high) )

N=$num_cores
(
for k in *.yml
do
    ((i=i%N)); ((i++==0)) && wait
    fname=${k%.*}
    export SOL_DIR="${cwd}/${run_dir}"
    export SOL_FILE="${SOL_DIR}/${fname}"
    FILE=${SOL_FILE}.h5
    if [ -f "$FILE" ]; then
        echo -e "${SOL_FILE} exists. Going to next config \n"
	    #mpirun -np 9 python pnp_par.py $SOL_FILE &
        continue
    else
        echo -e "Running config file : ${k}"
	    #cp ../${base_name}.yml ${SOL_FILE}.yml
        envsubst < ${k} > ${k}
        echo $SOL_FILE
        #echo -e "Running ${fname} \n"
        mpirun -np 9 python pnp_par.py $SOL_FILE && python get_current.py $SOL_FILE &
        #python get_current.py $SOL_FILE
        echo -e "Ended ${fname} \n"
    fi
done
)
#wait
echo "All done"
cd ../

end=`date +%s`
runtime=$((end-start))
echo -e "\n Runtime : ${runtime} \n"
