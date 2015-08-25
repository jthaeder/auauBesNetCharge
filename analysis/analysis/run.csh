#!/bin/csh

make clean
make 

./makeRefMultList.sh file.list 1

./analysis --energy=$1



