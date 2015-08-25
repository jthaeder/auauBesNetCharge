#!/bin/csh


make clean

#starver SL11d
make 

./makeRefMultList.sh file.list 1

./analysis --energy=11.5



