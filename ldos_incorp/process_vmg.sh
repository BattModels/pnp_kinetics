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
if [ -f "$out" ]; then
    echo "$out exists. Overwriting"
    echo "I Vapp Vdl" > $out
else
    echo "I Vapp Vdl" > $out
fi

Vapplist=()
N=$num_cores
(
for f in *.yml
do
    ((i=i%N)); ((i++==0)) && wait
    echo -e "\nAnalyzing ${f}"
    fname=${f%.*}
    count=`ls -1 ${fname}.h5 2>/dev/null | wc -l`
if [ $count != 0 ]; then
    echo "Output file exists ${fname}.h5, getting current..."
    Vapplist+=$Vapp
    python ../get_current.py $fname &
    echo -e "---------------------------------------------------------------------------------\n"
fi
done
)
#wait 
echo "All done"

cd ../