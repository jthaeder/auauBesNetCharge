#! /bin/tcsh
set nevents="-1"

set energy=$1

set particle=$2

set mode=1

if ( $mode == 0 ) then
    if ( $particle == 0 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/PiMinus_'$energy'.list.clean",    "PiMinus", '$energy')'
    else if ( $particle == 1 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/PiPlus_'$energy'.list.clean",     "PiPlus",  '$energy')'
    else if ( $particle == 2 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/AntiProton_'$energy'.list.clean", "Pbar",    '$energy')'
    else if ( $particle == 3 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/Proton_'$energy'.list.clean",     "Proton",  '$energy')'
    else if ( $particle == 4 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/KPlus_'$energy'.list.clean",      "KPlus",   '$energy')'
    else if ( $particle == 5 ) then
	root4star -b -q minimc_macro.cc'('$nevents',"lists/KMinus_'$energy'.list.clean",     "KMinus",  '$energy')'
    endif
endif

if ( $mode == 1 ) then
    if ( $particle == 0 ) then
	root4star -b -q make_efficiencies.C'("piplus",       '$energy')'
    else if ( $particle == 1 ) then
	root4star -b -q make_efficiencies.C'("piminus",      '$energy')'
    else if ( $particle == 2 ) then
	root4star -b -q make_efficiencies.C'("kaonplus",     '$energy')'
    else if ( $particle == 3 ) then
	root4star -b -q make_efficiencies.C'("kaonminus",    '$energy')'
    else if ( $particle == 4 ) then
	root4star -b -q make_efficiencies.C'("protonplus",   '$energy')'
    else if ( $particle == 5 ) then
	root4star -b -q make_efficiencies.C'("protonminus",  '$energy')'
    endif
endif
