#!/bin/bash

#for ii in 7.7 11.5 19.6 27 39 62.4 200 ; do
#    ./submit.sh $ii base
#done


for chargeSep in 1 2 ; do
#    for eta in 0.5 0.4 0.3 0.2 0.1 ; do 
    for eta in 0.45 0.35 0.25 0.15 0.05 ; do 
	./submit.sh 14.5 chargeSep_${chargeSep}_-${eta}_${eta}  ${chargeSep} -${eta} ${eta}
    done
done

