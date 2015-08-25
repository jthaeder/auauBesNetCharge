#!/bin/bash


inFile=file.list.15GeV
outFile=file.list.15GeV.clean

badrunFile=bad_runs_3sig.txt

tmpFile=file.list.15GeV.tmp

cp $inFile $tmpFile
count=0
while read -r line ; do
    let count=count+1
    echo "Process [$count] : $line"
    grep -v $line $tmpFile > $outFile 
    cp $outFile $tmpFile
done < <(cat $badrunFile | sort | uniq )

rm $tmpFile 

