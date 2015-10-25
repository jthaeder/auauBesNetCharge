#!/bin/bash

energies="200"

mkdir -p log embeddingTrees efficiency efficiencyStudy

for energy in $energies ; do 
    for particle in 0 1 2 3 4 5  ; do 
	./run.csh $energy $particle 2>&1 log/job_${energy}_${particle}.err
    done
done

    
