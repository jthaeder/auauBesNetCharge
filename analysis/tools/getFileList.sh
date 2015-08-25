#!/bin/bash


module use /common/star/pkg/Modules
module load python pymongo star-dm-scripts



for energy in 7.7GeV 11GeV 39GeV 62GeV 200GeV ; do 

    list=energy_fileLists/filelist_${energy}.list
    if [ "$energy" = "200GeV" ] ; then 
	list=energy_fileLists/filelist_${energy}_Run10.list
    fi

    if [ -f $list ] ; then 
	rm $list 
    fi

    starquery '{"runyear" : "Run10",  "energy" : "'$energy'", "stream" : {"$regex" : "st_physics.*"}}' | sort | uniq | sed -e 's|^|root://pstarxrdr1/|' > $list
    
    cat $list | cut -d'/' -f 9 | sort | uniq
    cat $list | cut -d'/' -f 10 | sort | uniq
    cat $list | cut -d'/' -f 11 | sort | uniq
done
    

for energy in 19GeV 27GeV 200GeV ; do 
    list=energy_fileLists/filelist_${energy}.list
    if [ "$energy" = "200GeV" ] ; then 
	list=energy_fileLists/filelist_${energy}_Run11.list
    fi
    if [ -f $list ] ; then 
	rm $list 
    fi

    starquery '{"runyear" : "Run11",  "energy" : "'$energy'", "stream" : {"$regex" : "st_physics.*"}}' | sort | uniq | sed -e 's|^|root://pstarxrdr1/|' > $list
    
    cat $list | cut -d'/' -f 9 | sort | uniq
    cat $list | cut -d'/' -f 10 | sort | uniq
    cat $list | cut -d'/' -f 11 | sort | uniq
done