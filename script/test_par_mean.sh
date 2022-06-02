#!/bin/bash
# This bash script lauches the input program $1 and tests it with increasing parallelism degree in a loop up to 'par'.
# Each iteration is an average of 'iter' number of independent runs.

if [ $# -ne 4 ]; then
    echo "./bin/test [matrix_dim_value] [method] [max_nw] [iter]"
	exit 1
fi

n=$1		 	# matrix dimension nxn
m=$2			# method to use [seq, th, ff]
max_nw=$3		# maximum parallelism degree nw
iter=$4			# number of iterations to reduce variance



echo "nw|   average"

# loop to increase the parallelism degree
for i in $(seq 1 1 $max_nw)
do
	# loop for multiple independent iterations to reduce variance
	for j in  $(seq 1 1 $iter)
	do
		./bin/test -n $n -m $m -w $i -l -2 -r 2 -t 10e-6 -s 20; 
	done | grep jacobi | awk -v iter=$iter -v par=$i '{sum += $4} END {print sprintf("%s | \t%.3f", par, sum/iter)}'
done