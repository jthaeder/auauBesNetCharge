#!/bin/bash

particles="KPlus KMinus AntiProton Proton PiPlus PiMinus"

outputdir=`pwd`/lists

for energy in 11 14 19 ; do 
    for p in ${particles} ; do 
	
	outfile=${outputdir}/${p}_${energy}.list
	
	outfileClean=${outputdir}/${p}_${energy}.list.clean
	
	if [ -f $outfileClean ] ; then
	    rm $outfileClean
	fi
	
	touch   $outfileClean
       	
	while read -r line ; do 
	    if [ -s $line ] ; then
		grep $line $outputdir/badFiles.list > /dev/null
		if [ $? -ne 0 ] ; then
		    echo $line >> $outfileClean
		fi
	    fi
	done < <(cat $outfile)
	
	lines=`cat $outfile | wc -l`
	linesClean=`cat $outfileClean | wc -l`
	
	echo $outfile $lines $linesClean
    done
done
