#!/bin/bash


# Usage : makeRefMultList.sh <fileList> (<mode>) 

# mode:
#   0 - make refMultList for filelist (default)
#   1 - clean file list
mode=0


fileList=$1
if [ $# -eq 0 ] ; then 
    exit
elif [ $# -eq 2 ] ; then 
    mode=$2
fi

fileListHasRefMult=${fileList}.hasRefmult
if [ -f ${fileListHasRefMult} ] ; then
    rm ${fileListHasRefMult}
fi

fileListHasNoRefMult=${fileList}.hasNoRefmult
if [ -f ${fileListHasNoRefMult} ] ; then
    rm ${fileListHasNoRefMult}
fi

fileListRefMult=${fileList}.refMult.tmp
if [ -f ${fileListRefMult} ] ; then
    rm ${fileListRefMult}
fi

while read -r file ; do 
    energy=`echo $file | cut -d'/' -f 9 | awk -F'GeV' '{ print $1 }'`
    day=`echo $file | cut -d'/' -f 12`
    run=`echo $file | cut -d'/' -f 13`
    name=`echo $file | cut -d'/' -f 14 | awk -F'.' '{ print $1 }'`

    if [ "$energy" = "7.7" ] ; then 
	energy=7
    fi

    refMultFile=refMultExtract/refMult_${energy}/${day}/${run}/${name}.refMult.txt

    if [ -f ${refMultFile} ] ; then 
	echo ${file} >> ${fileListHasRefMult}
	if [ $mode -eq 1 ] ; then
	    echo ${refMultFile} >> ${fileListRefMult}
	fi
    else
	echo ${file} >> ${fileListHasNoRefMult}
    fi
done < <(cat $fileList)


if [ $mode -eq 1 ] ; then

    listRefMult=${fileList}.refMult
    if [ -f ${listRefMult} ] ; then
	rm ${listRefMult}
    fi
    
    name=`cat $fileListRefMult | head -n 1`
    cat $name | head -n 1 >> ${listRefMult}
    
    while read -r file ; do 
	cat $file | grep -v event >> ${listRefMult}
    done < <(cat $fileListRefMult)

    rm -f $fileListRefMult
fi
