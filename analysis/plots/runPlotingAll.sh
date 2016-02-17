#!/bin/bash

rm -f STAR_QM2015_Preliminary.root STAR_Preliminary.root
rm -f *.so *.d include/*.so include/*.d

./runPlotingEnergy.sh
./runPlotingEta.sh

root -l -b -q makeBESIIerror.C++
root -l -b -q plotCumulants.C++
#root -l -b -q plotSetNetProtonCrossCheck.C
#root -l -b -q plotSetSys.C


# -- Create QM 2015 data files
root -l -b -q print_STAR_QM2015_Preliminary.C++ > STAR_QM2015_Preliminary.txt

rm -f STAR_QM2015_Preliminary.zip

zip STAR_QM2015_Preliminary.zip STAR_QM2015_Preliminary.root STAR_QM2015_Preliminary.txt print_STAR_QM2015_Preliminary.C

rm -f *.so *.d include/*.so include/*.d
