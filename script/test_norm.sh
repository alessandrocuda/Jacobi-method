#!/bin/bash
# Each iteration is an average of 'iter' number of independent runs.

if [ $# -ne 1 ]; then
            echo "./bin/test [iter]"
                exit 1
fi

iter=$1         # number of iterations to reduce variance



echo "matrix dim |   average     | norm % | norm ideal %"

# loop to increase the parallelism degree
for i in 4 16 64 128 500 1000 5000 10000 20000 30000 50000
do
        # loop for multiple independent iterations to reduce variance
        for j in  $(seq 1 1 $iter)
        do
        ./bin/test -n $i -m seq;
done | grep norm | awk -v iter=$iter -v par=$i '{sum += $3} {sum2 += $6} END {print sprintf("%s | \t%.3f | \t %.3f | \t %.3f | \t %.3f", par, sum/iter, sum2/iter, (100/sum)*sum2, (100/(par*par +par))*par)}'
done