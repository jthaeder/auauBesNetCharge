#!/bin/sh

get_file_list.pl -keys 'runnumber' -cond 'production=P14ii,trgsetupname=production_15GeV_2014,filename~st_physics,' -limit 0 -distinct | sort > run.list

