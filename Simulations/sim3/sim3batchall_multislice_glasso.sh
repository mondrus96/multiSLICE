#!/bin/bash
iters=(1 2 3 4)
ns=(25 50 100 200) #sample size and p size
Sest="glasso"

# Loop through all combos
for iter in "${iters[@]}"; do
    for n in "${ns[@]}"; do
        # sbatch --job-name="${Sest}_n${n}_l${l}_iter${iter}" --time=0-11:59 sim1batch.sh $iter $n $Sest
        Rscript runsim3.R $iter $n $Sest
        sleep 1
    done
done