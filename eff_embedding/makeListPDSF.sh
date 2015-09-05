#!/bin/bash


particles="KPlus KMinus AntiProton Proton PiPlus PiMinus"


outputdir=`pwd`/lists

#
# make Lists
#


for energy in 7 11 14 19 27 39 62 200; do
    if [ $energy -eq 7 ] ; then
	basedirString=AuAu7_production
    elif [ $energy -eq 11 ] ; then
	basedirString=AuAu11_production
    elif [ $energy -eq 14 ] ; then
	basedirString=production_15GeV_2014
    elif [ $energy -eq 19 ] ; then
	basedirString=AuAu19_production
    elif [ $energy -eq 27 ] ; then
	basedirString=AuAu27_production_2011
    elif [ $energy -eq 39 ] ; then
	basedirString=AuAu39_production
    elif [ $energy -eq 62 ] ; then
	basedirString=AuAu62_production
    elif [ $energy -eq 200 ] ; then
	basedirString=2010ProductionMinBias
        #AuAu200_production
    fi

    for p in ${particles} ; do 

	outfile=${outputdir}/${p}_${energy}.list.pdsf

	if [ -f $outfile ] ; then
	    rm $outfile
	fi
	
	touch $outfile
	touch ${outfile}.tmp

	basedir=/projecta/projectdirs/starprod/embedding/${basedirString}
	echo $basedir
	if [ ! -d $basedir ] ; then
	    continue;
	fi
	
	pushd $basedir > /dev/null
	
	for folder in `ls | grep -i $p` ; do 
	    find ${basedir}/${folder}  -iname "*minimc.root" >> ${outfile}.tmp
	done
	
	if [ "$p" = "AntiProton" ] ; then
	    for folder in `ls | grep Pbar` ; do 
		find ${basedir}/${folder}  -iname "*minimc.root" >> ${outfile}.tmp
	    done
	fi
	popd > /dev/null
		
	cat  ${outfile}.tmp| sort | uniq >  ${outfile}
	rm ${outfile}.tmp
	
	nFiles=`cat ${outfile} | wc -l`
	if [ $nFiles -eq 0 ] ; then 
	    rm ${outfile}
	fi
    done
done


