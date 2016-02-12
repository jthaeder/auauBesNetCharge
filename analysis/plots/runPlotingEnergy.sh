#!/bin/bash

rm -f STAR_QM2015_Preliminary.root STAR_Preliminary.root
rm -f *.so *.d include/*.so include/*.d

root -l -b -q calcSys.C++

root -l -b -q plotEnergyCharge.C++
root -l -b -q plotEnergyProton.C++
root -l -b -q plotEnergyKaon.C++
root -l -b -q plotEnergyProtonOverview.C++

root -l -b -q plotRatioSummary.C++

rm -f *.so *.d include/*.so include/*.d

