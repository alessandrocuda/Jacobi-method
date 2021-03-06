#!/bin/bash
echo "matrix dim |   T_jacobi_seq |     T_error_tot     |       p        | Ideal p      "

# loop to increase the parallelism degree
for i in 4 16 64 128 500 1000 5000 10000 20000 30000 50000
do
        # loop for multiple independent iterations to reduce variance
        for j in  $(seq 1 1 10)
        do
        ./bin/test -n $i -m seq -s 20 -l -2 -r 2 -t 10e-6 -v 1;
done | grep norm | awk -v iter=10 -v par=$i '{sum += $3} {sum2 += $6} END {print sprintf("%s | \t%.3f | \t %.3f | \t %.3f | \t %.3f", par, sum/iter, sum2/iter, (100/sum)*sum2, (100/(par*par +par))*par)}'
done