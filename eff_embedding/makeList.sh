#!/bin/bash


particles="KPlus KMinus AntiProton Proton PiPlus PiMinus"


outputdir=`pwd`/lists

#
# make Lists
#

disks="data18 data19 data20 data21 data22"


for energy in 7 11 14 19 27 39 62 ; do
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
	basedirString=AuAu200_production
    fi

    for p in ${particles} ; do 

	outfile=${outputdir}/${p}_${energy}.list

	if [ -f $outfile ] ; then
	    rm $outfile
	fi
	
	touch $outfile
	touch ${outfile}.tmp

	for ii in $disks ; do
	    basedir=/star/${ii}/embedding/${basedirString}
		
	    if [ ! -d $basedir ] ; then
		continue;
	    fi

	    pushd $basedir > /dev/null
	    
	    for folder in `ls | grep -i $p` ; do 
		find ${basedir}/${folder}  -name "*minimc.root" >> ${outfile}.tmp
	    done

	    if [ "$p" = "AntiProton" ] ; then
		for folder in `ls | grep Pbar` ; do 
		    find ${basedir}/${folder}  -name "*minimc.root" >> ${outfile}.tmp
		done
	    fi
	    popd > /dev/null
	done
	
	cat  ${outfile}.tmp| sort | uniq >  ${outfile}
	rm ${outfile}.tmp
    done
done

exit

#
# 11 GEV
#

for p in ${particles} ; do 

    outfile=${outputdir}/${p}_11.list

    if [ -f $outfile ] ; then
	rm $outfile
    fi

    touch $outfile
    touch ${outfile}.tmp


    for ii in data21 ; do
	
	basedir=/star/${ii}/embedding/AuAu11_production

	pushd $basedir > /dev/null
	
	for folder in `ls | grep -i $p` ; do 
	    find ${basedir}/${folder}  -name "*minimc.root" >> ${outfile}.tmp
	done
	popd > /dev/null
    done
    
    cat  ${outfile}.tmp| sort | uniq >  ${outfile}
    rm ${outfile}.tmp
done


#
# 19 GEV
#

for p in ${particles} ; do 

    outfile=${outputdir}/${p}_19.list

    if [ -f $outfile ] ; then
	rm $outfile
    fi

    touch $outfile
    touch ${outfile}.tmp


    for ii in data19 ; do
	
	basedir=/star/${ii}/embedding/AuAu19_production

	pushd $basedir > /dev/null
	
	for folder in `ls | grep -i $p` ; do 
	    find ${basedir}/${folder}  -name "*minimc.root" >> ${outfile}.tmp
	done
	popd > /dev/null
    done
    
    cat  ${outfile}.tmp| sort | uniq >  ${outfile}
    rm ${outfile}.tmp
done


