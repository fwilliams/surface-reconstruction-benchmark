#!/bin/bash

# format: ./poisson_recon.sh in_file out_file remaining params
in_file=$1
out_file=$2

poisson_args="./recon/poisson/PoissonRecon --in $in_file --out $out_file"

arg_ind=0
for arg in $@
do
	let arg_ind=$arg_ind+1
	if [[ $arg_ind -lt 3 ]]
	then
		continue
	fi
	poisson_args="$poisson_args $arg"
done

std_out=${out_file%.ply}.out
std_err=${out_file%.ply}.err

echo "poisson reconstruction..."
$poisson_args > $std_out 2> $std_err
echo "... poisson reconstruction done"
