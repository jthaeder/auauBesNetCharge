#!/bin/bash

./clean.sh

root -l -b -q fitEfficiencyPt.C++
root -l -b -q fitEfficiencyEta.C++
root -l -b -q makeEfficiencyCharged.C++
#root -l -b -q makeEfficiencyStudy.C++

rm -f *.d *.so
