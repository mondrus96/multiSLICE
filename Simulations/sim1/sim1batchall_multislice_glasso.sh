#!/bin/bash
iters=(1 2 3 4)
ns=(100 200 300 400 500) #sample size
ps=(2 3 4 5) #layers
Sest="glasso"

# Loop through all combos
for iter in "${iters[@]}"; do
    for n in "${ns[@]}"; do
        for l in "${ls[@]}"; do
            sbatch --job-name="${Sest}_n${n}_l${l}_iter${iter}" --time=0-11:59 sim1batch.sh $iter $n $l $Sest
            sleep 1
        done
    done
done