#!/bin/bash

energies="11 14 19 27 39 62 7"
# 200

for energy in $energies ; do 
    echo "Arguments      = ${energy}" > submit_condor.con
    cat submit_condor.template >> submit_condor.con

    condor_submit submit_condor.con
    rm  submit_condor.con
done

    
