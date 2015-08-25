#!/bin/csh

set inFile=$1

rm -rf filelist datalist

root -l -b -q Split.C'("'${inFile}'")'
