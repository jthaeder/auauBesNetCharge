#!/bin/bash

##  ./submit.sh <energy> <suffix>

energy=$1
suffix=$2

folder=jobs

if [ $# -eq 2 ] ; then
    folder=jobs/${folder}_${energy}_${suffix}
elif [ $# -eq 1 ] ; then
    folder=jobs/${folder}_${energy}
else
    echo  "   ./submit.sh <energy> <suffix>"
    exit 0
fi

if [ "$energy" == "7.7" ] ; then
    fileList=filelist_7.7GeV.list.hasRefmult
elif [ "$energy" == "11.5" ] ; then
    fileList=filelist_11GeV.list.hasRefmult
elif [ "$energy" == "14.5" ] ; then
    fileList=filelist_14GeV.list.clean
elif [ "$energy" == "19.6" ] ; then
    fileList=filelist_19GeV.list.hasRefmult
elif [ "$energy" == "27" ] ; then
    fileList=filelist_27GeV.list.hasRefmult
elif [ "$energy" == "39" ] ; then
    fileList=filelist_39GeV.list.hasRefmult
elif [ "$energy" == "62.4" ] ; then
    fileList=filelist_62GeV.list.hasRefmult
else
    echo  "   ./submit.sh <energy> <suffix>"
    exit 0
fi

./runSplit.csh energy_fileLists/${fileList}

jobIdx=0

for file in `cat datalist` ; do 
    rm -rf $folder/$jobIdx
    mkdir -p $folder

    cp -a analysis $folder/$jobIdx

    pushd $folder/$jobIdx > /dev/null
    rm file.list*
    cp ../../../$file ./file.list

    qsub -l h_vmem=3G -l projectio=1,scratchfree=500,h_cpu=24:00:00 -o ./job.out -e ./job.err run.csh --energy=${energy}
    popd > /dev/null

    let "jobIdx+=1"
done
