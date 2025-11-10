#!/bin/bash
iters=(1 2 3 4)
ns=(100 200 300 400 500) #sample size
plats=(2 3 4 5) #plat
Sest="glasso"

# Loop through all combos
for iter in "${iters[@]}"; do
    for n in "${ns[@]}"; do
        for plat in "${plats[@]}"; do
            sbatch --job-name="${Sest}_n${n}_plat${plat}_iter${iter}" --time=0-11:59 sim2batch.sh $iter $n $plat $Sest
            sleep 1
        done
    done
done