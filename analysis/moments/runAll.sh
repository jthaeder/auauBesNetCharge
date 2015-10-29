#!/bin/bash

mode=2

# -- sys studies
#versions="jobs_14.5_base jobs_14.5_sys_0_0 jobs_14.5_sys_0_1 jobs_14.5_sys_0_2 jobs_14.5_sys_0_3 jobs_14.5_sys_1_0 jobs_14.5_sys_1_1 jobs_14.5_sys_1_2 jobs_14.5_sys_1_3"

# ------------------------------------------------------------------

#versions="jobs_14.5_base_plus5"
#versions="jobs_14.5_base_minus5"
#versions="jobs_14.5_base_plus5minus5"
#versions="jobs_14.5_base_plus2minus2"
#versions="jobs_14.5_base_minus2plus2"
#versions="jobs_11.5_test"
#versions="jobs_19.6_test"
# ------------------------------------------------------------------

# -- all delta eta studies
versions="2015-05-20_0.1 2015-05-20_0.2 2015-05-20_0.3 2015-05-20_0.4 2015-05-20_0.5 2015-06-03_delta_0.3_0 2015-06-03_delta_0.3_1 2015-06-03_delta_0.3_2 2015-06-03_delta_0.3_3 2015-06-03_delta_0.3_4 2015-06-03_delta_0.3_5 2015-06-03_delta_0.3_6 2015-06-03_delta_0.3_7 2015-06-03_delta_0.3_8 2015-06-03_delta_0.5_0 2015-06-03_delta_0.5_1 2015-06-03_delta_0.5_2 2015-06-03_delta_0.5_3 2015-06-03_delta_0.5_4 2015-06-03_delta_0.5_5 2015-06-20_delta_0.5_6 2015-06-20_delta_0.7_0 2015-06-20_delta_0.9_0 2015-06-21_delta_0.1_0 2015-06-21_delta_1.0_mult_0.5"

# -- mid rapidty (delta eta [0.1, 1.0])
#versions="2015-06-21_delta_0.1_0 2015-05-20_0.1 2015-06-03_delta_0.3_8 2015-05-20_0.2 2015-06-20_delta_0.5_6 2015-05-20_0.3 2015-06-20_delta_0.7_0 2015-05-20_0.4 2015-06-20_delta_0.9_0 2015-05-20_0.5"

# -- mid rapidty (delta eta [0.1, 1.0]) - with charge separation
#versions="jobs_14.5_chargeSep_1_-0.05_0.05 jobs_14.5_chargeSep_1_-0.1_0.1 jobs_14.5_chargeSep_1_-0.15_0.15 jobs_14.5_chargeSep_1_-0.2_0.2 jobs_14.5_chargeSep_1_-0.25_0.25 jobs_14.5_chargeSep_1_-0.3_0.3 jobs_14.5_chargeSep_1_-0.35_0.35 jobs_14.5_chargeSep_1_-0.4_0.4 jobs_14.5_chargeSep_1_-0.45_0.45 jobs_14.5_chargeSep_1_-0.5_0.5 jobs_14.5_chargeSep_2_-0.05_0.05 jobs_14.5_chargeSep_2_-0.1_0.1 jobs_14.5_chargeSep_2_-0.15_0.15 jobs_14.5_chargeSep_2_-0.2_0.2 jobs_14.5_chargeSep_2_-0.25_0.25 jobs_14.5_chargeSep_2_-0.3_0.3 jobs_14.5_chargeSep_2_-0.35_0.35 jobs_14.5_chargeSep_2_-0.4_0.4 jobs_14.5_chargeSep_2_-0.45_0.45 jobs_14.5_chargeSep_2_-0.5_0.5"

# ------------------------------------------------------------------

dataSets="effuncorr twoeff_11"
dataSets="twoeff_11"
dataSets="twoeff_11_pos twoeff_11_neg"


name=Sum_Moments_hist_AuAu14.5GeV_charge.root
#name=Sum_Moments_hist_AuAu11.5GeV_charge.root
#name=Sum_Moments_hist_AuAu19.6GeV_charge.root

for version in $versions ; do

    mkdir -p output/${version}/
    
    # -- make fact4
    # -------------------------
    if [[ $mode -eq 0 || $mode -eq 2 ]] ; then
	rm -f output/${version}/Fact4.root
	./runMakeFact4.csh ${version} ${name}
    fi
    
    for dataSet in $dataSets ; do 

	mkdir -p output/${version}/${dataSet}
	pushd output/${version}/${dataSet} > /dev/null
	if [ ! -h Fact4.root ] ; then
	    ln -s ../Fact4.root
	fi
	popd > /dev/null
	
	# -- make fact
	# -------------------------
	if [[ $mode -eq 0 || $mode -eq 2 ]] ; then
	    pushd Convert2FF_${dataSet} > /dev/null
	    rm -f Fact.root analysis.o analysis

	    ln -sf ../output/${version}/${dataSet}/Fact4.root

	    ./run.csh 2>&1 | tee ../output/${version}/${dataSet}/convert2Fact.log
	    cp Fact.root ../output/${version}/${dataSet}/

	    popd > /dev/null
	fi
	
	if [[ $mode -eq 1 || $mode -eq 2 ]] ; then
	    pushd moments_ana_${dataSet}  > /dev/null
	    rm -f Moments_*.root analysis.o analysis

	    ln -sf ../output/${version}/${dataSet}/Fact4.root
	    ln -sf ../output/${version}/${dataSet}/Fact.root

	    ./run.csh 2>&1 | tee ../output/${version}/${dataSet}/moments.log
	    cp Moments_*.root ../output/${version}/${dataSet}/

	    popd > /dev/null
	fi
    done
done
