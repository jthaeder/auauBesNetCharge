#!/bin/bash

energies="7 11 14 19 27 39 62 200"

mkdir -p log embeddingTrees efficiency efficiencyStudy

for energy in $energies ; do 
    for particle in 0 1 2 3 4 5  ; do 
	./run.csh ${energy} ${particle}
    done
done

    
