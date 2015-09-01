#!/bin/bash

energies="7 11 14 19 27 39 62"
# 200

mkdir -p log

for energy in $energies ; do 
    for particle in 0 1 2 3 4 5  ; do 
	echo "Arguments      = ${energy} ${particle}" > submit_condor.con
	echo "Log            = log/job_${energy}_${particle}.log" >> submit_condor.con
	echo "Output         = log/job_${energy}_${particle}.out" >> submit_condor.con
	echo "Error          = log/job_${energy}_${particle}.err" >> submit_condor.con
	
	cat submit_condor.template >> submit_condor.con

	condor_submit submit_condor.con
	rm  submit_condor.con
    done
done

    
