#!/bin/csh

make clean
make 

set energy=$1

if ( $energy != 14.5 ) then
    ./makeRefMultList.sh file.list 1
endif

set chargeSeparation=$2
set etaMin=$3
set etaMax=$4

./analysis --energy=$energy --chargeSeparation=$chargeSeparation --etaMin=$etaMin --etaMax=$etaMax
