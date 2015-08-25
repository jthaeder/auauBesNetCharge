#!/bin/csh

make clean
make 

set energy=$1

if ( $energy != 14.5 ) then
    ./makeRefMultList.sh file.list 1
endif

./analysis --energy=$energy



