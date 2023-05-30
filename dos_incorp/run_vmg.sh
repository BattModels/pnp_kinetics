#!/bin/bash

num_cores=3
out=iv_data.txt
dir=$1
if [ -z "$1" ]; then
    echo "No argument supplied, using dir=Sol/"
    dir=Sol
    if [ ! -d $dir ]; then
    echo "Run directory ${dir} does not exist. Making new dir Sol/"
    mkdir -p $dir
    ./prep_vmg.sh $dir
    fi
else
    dir=$1
fi

cd $dir
FILE=iv_data.txt
if [ -f "$FILE" ]; then
    echo "$FILE exists. Appending file"
else
    echo "I Vapp Vdl" > $out
fi

for f in *.yml
do
    echo -e "\nAnalyzing ${f}"
    fname=${f%.*}
    count=`ls -1 ${fname}.h5 2>/dev/null | wc -l`
if [ $count != 0 ]; then 
    echo "Output file already exists ${fname}.h5, getting current..."
    python ../get_current.py $fname
    echo -e "---------------------------------------------------------------------------------\n"
continue
else
    echo -e "Running ${f} \n"
    mpirun -np $num_cores python ../pnp_par.py $fname
    python ../get_current.py $fname 
    echo "Ended ${f}"
    echo -e "---------------------------------------------------------------------------------\n"
fi
done
cd ../
