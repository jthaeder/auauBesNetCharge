#!/bin/bash

./runPlotingEnergy.sh
./runPlotingEta.sh

root -l -b -q makeBESIIerror.C++
root -l -b -q plotCumulants.C++
#root -l -b -q plotSetNetProtonCrossCheck.C
#root -l -b -q plotSetSys.C

root -l -b -q print_STAR_QM2015_Preliminary.C++ > STAR_QM2015_Preliminary.txt

rm -f *.so *.d include/*.so include/*.d
