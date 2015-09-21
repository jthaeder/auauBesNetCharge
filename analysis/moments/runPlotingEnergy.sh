#!/bin/bash

rm -f *.so *.d include/*.so include/*.d

root -l -b -q calcSys.C++

root -l -b -q plotEnergyCharge.C++
root -l -b -q plotEnergyProton.C++
root -l -b -q plotEnergyKaon.C++

rm -f *.so *.d include/*.so include/*.d

