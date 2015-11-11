#!/bin/csh

make clean
make 

set analysis=$1
set energy=$2

if ( $energy != 14.5 ) then
    ./makeRefMultList.sh file.list 1
endif

set chargeSeparation=$3
set etaMin=$4
set etaMax=$5

./analysis --analysis=$analysis --energy=$energy 
# --chargeSeparation=$chargeSeparation --etaMin=$etaMin --etaMax=$etaMax
