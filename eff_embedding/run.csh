#! /bin/tcsh
set nevents="-1"

set energy=$1

set particle=$2

set mode=2

if ( $mode == 0 ) then
    if ( $particle == 0 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/PiMinus_'$energy'.list.pdsf",    "PiMinus", '$energy')'
    else if ( $particle == 1 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/PiPlus_'$energy'.list.pdsf",     "PiPlus",  '$energy')'
    else if ( $particle == 2 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/AntiProton_'$energy'.list.pdsf", "Pbar",    '$energy')'
    else if ( $particle == 3 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/Proton_'$energy'.list.pdsf",     "Proton",  '$energy')'
    else if ( $particle == 4 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/KPlus_'$energy'.list.pdsf",      "KPlus",   '$energy')'
    else if ( $particle == 5 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/KMinus_'$energy'.list.pdsf",     "KMinus",  '$energy')'
    endif
endif

if ( $mode == 1 ) then
    if ( $particle == 0 ) then
	root4star -b -q make_efficiencies.C+'("piplus",       '$energy')'
    else if ( $particle == 1 ) then
	root4star -b -q make_efficiencies.C+'("piminus",      '$energy')'
    else if ( $particle == 2 ) then
	root4star -b -q make_efficiencies.C+'("kaonplus",     '$energy')'
    else if ( $particle == 3 ) then
	root4star -b -q make_efficiencies.C+'("kaonminus",    '$energy')'
    else if ( $particle == 4 ) then
	root4star -b -q make_efficiencies.C+'("protonplus",   '$energy')'
    else if ( $particle == 5 ) then
	root4star -b -q make_efficiencies.C+'("protonminus",  '$energy')'
    endif
endif

if ( $mode == 2 ) then
    if ( $particle == 0 ) then
	root4star -b -q study_efficiencies.C+'("piplus",       '$energy')'
    else if ( $particle == 1 ) then
	root4star -b -q study_efficiencies.C+'("piminus",      '$energy')'
    else if ( $particle == 2 ) then
	root4star -b -q study_efficiencies.C+'("kaonplus",     '$energy')'
    else if ( $particle == 3 ) then
	root4star -b -q study_efficiencies.C+'("kaonminus",    '$energy')'
    else if ( $particle == 4 ) then
	root4star -b -q study_efficiencies.C+'("protonplus",   '$energy')'
    else if ( $particle == 5 ) then
	root4star -b -q study_efficiencies.C+'("protonminus",  '$energy')'
    endif
endif
