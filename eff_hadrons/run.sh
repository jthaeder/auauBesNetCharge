#!/bin/bash

./clean.sh

root -b -l -q fitEfficiencyPt.C++
root -b -l -q fitEfficiencyEta.C++
root -l -b -q makeEfficiencyCharged.C++

rm -f *.d *.so
