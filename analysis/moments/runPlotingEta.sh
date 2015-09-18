#!/bin/bash

rm -f *.so *.d

root -l -b -q plotEta.C++
root -l -b -q plotVsDeltaEta.C++
root -l -b -q plotVsEta.C++

root -l -b -q plotEta_nice.C++
root -l -b -q plotVsDeltaEta_nice.C++
root -l -b -q plotVsEta_nice.C++

rm -f *.so *.d
