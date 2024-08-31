#!/bin/bash

start=`date +%s`

mpirun -np "$1" python ../pnp_par.py "$2" && python update_current.py "$2"

end=`date +%s`
runtime=$((end-start))
echo -e "\n Runtime : ${runtime} \n"