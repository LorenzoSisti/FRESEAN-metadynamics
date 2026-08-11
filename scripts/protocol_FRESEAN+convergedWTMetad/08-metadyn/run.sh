#!/bin/bash

nReplicas=20
for (( i=1; i<${nReplicas}; i++ ))
do
cp -r single_metad metadyn_$i
cd metadyn_$i
echo "Run replica $i"
sbatch --job-name=METAD_REPLICA_${i}.run --export=replica=${i} startme.sh
cd ..
done
