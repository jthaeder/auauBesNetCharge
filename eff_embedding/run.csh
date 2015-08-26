#! /bin/tcsh
set nevents="-1"

set energy=$1

root4star -b -q minimc_macro.cc'('$nevents',"lists/PiMinus_'$energy'.list.clean",    "PiMinus", '$energy')'
root4star -b -q minimc_macro.cc'('$nevents',"lists/PiPlus_'$energy'.list.clean",     "PiPlus",  '$energy')'

root4star -b -q minimc_macro.cc'('$nevents',"lists/AntiProton_'$energy'.list.clean", "Pbar",    '$energy')'
root4star -b -q minimc_macro.cc'('$nevents',"lists/Proton_'$energy'.list.clean",     "Proton",  '$energy')'

root4star -b -q minimc_macro.cc'('$nevents',"lists/KPlus_'$energy'.list.clean",      "KPlus",   '$energy')'
root4star -b -q minimc_macro.cc'('$nevents',"lists/KMinus_'$energy'.list.clean",     "KMinus",  '$energy')'

exit

root4star -b -q make_efficiencies.C'("piplus",       '$energy')'
root4star -b -q make_efficiencies.C'("piminus",      '$energy')'

root4star -b -q make_efficiencies.C'("kaonplus",     '$energy')'
root4star -b -q make_efficiencies.C'("kaonminus",    '$energy')'

root4star -b -q make_efficiencies.C'("protonplus",   '$energy')'
root4star -b -q make_efficiencies.C'("protonminus",  '$energy')'
