#!/bin/bash

particles="KPlus KMinus AntiProton Proton PiPlus PiMinus"

outputdir=`pwd`/lists


for energy in 7 11 14 19 27 39 62 ; do 
    cat log/job_${energy}*.err | grep " file " | grep "/star/" | awk -F '/star/' '{ print "/star/"$2 }' |  awk -F '.root' '{ print $1".root" }' > $outputdir/badFiles.list

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

    rm -f $outputdir/badFiles.list
done
