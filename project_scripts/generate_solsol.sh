#../aiida_create_solutesupercell_structures.py -a 4.04 -me Al -spsh 4,4,4 -fse Mg,Cu,Zn -sse Mg,Cu,Zn -sg SoluteSolute_structures -mxi 8
#../aiida_create_solutesupercell_structures.py -a 4.04 -me Al -spsh 4,4,4 -fse Mg,Si,Cu,Vac -sse Mg,Si,Cu,Vac -sg SoluteSolute_structures -mxi 8 
#../aiida_create_solutesupercell_structures.py -a 3.62 -me Cu -spsh 3,3,3 -fse Al,Vac -sse Al,Vac  -sg SoluteSolute_structures -mxi 4

#ELEMENTS=Li,Na,Ca,Sc,Ti,Cr,Mn,Co,Ni,Ge,Sr,Y,Zr,Pd,Ag,In,Sn,Sb,La,Pt,Au,Ce,Nd,Gd,Dy
#../aiida_create_solutesupercell_structures.py -a 4.04 -me Al -spsh 3,3,3 -fse $ELEMENTS -sse Cu -sg ManySolInAl_structures --single_solute_only

#../aiida_create_solutesupercell_structures.py -a 4.04 -me Al -spsh 4,4,4 -fse Zn -sse Zn,Vac,Mg,Cu -sg SoluteSolute_structures -mxi 4

#For Dilute Zn example
../aiida_create_solutesupercell_structures.py -a 4.04 -me Al -spsh 3,3,3 -fse Zn -sse Zn -sg DiluteZnInAl_structures --single_solute_only
