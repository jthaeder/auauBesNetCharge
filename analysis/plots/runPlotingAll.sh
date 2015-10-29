#!/bin/bash

./runPlotingEnergy.sh
./runPlotingEta.sh

root -l -b -q print_STAR_QM2015_Preliminary.C++ > STAR_QM2015_Preliminary.txt

rm -f *.so *.d include/*.so include/*.d
