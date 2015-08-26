#!/bin/bash

# Usage $0 <minIdx> <maxIdx>


for name in `ls jobs` ; do 
    outFolder=data/${name}
    
    if [ -d ${outFolder} ] ; then 
	continue
    fi

    mkdir -p $outFolder

    mergeList=${outFolder}/mergeRootFiles.list
    if [ -f $mergeList ] ; then
	rm -f $mergeList
    fi

    find jobs/$name -name "Moments_hist*.root" -type f  >> ${mergeList}
    
    fileName=`cat ${mergeList} | head -n 1`
    fileName=`basename $fileName`

    hadd -f ${outFolder}/Sum_${fileName} `cat ${mergeList}`
done




