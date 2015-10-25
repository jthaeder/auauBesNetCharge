#!/bin/csh


make clean

#starver SL11d
make 

#./makeRefMultList.sh file.list 1

./analysis --energy=14.5 --chargeSeparation=1 --etaMin=-0.25 --etaMax=0.25


