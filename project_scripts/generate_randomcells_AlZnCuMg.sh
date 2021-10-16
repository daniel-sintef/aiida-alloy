MATRIX_ELEMENTS=Al,Zn,Cu,Mg,Vac
MATRIX_LATTICES=4.07,3.02,3.65,4.58,-1
SG=Random-108-AlZnCuMg_structures
#DR=--dryrun
DR=
#NOSTORE=--nostore
NOSTORE=

CONCENTRATIONS=0.7,0.1,0.1,0.1,0.01
../aiida_create_randomsupercell_structures_v2.py -me $MATRIX_ELEMENTS  -c $CONCENTRATIONS  -ns 1000 -spsh 3,3,3 -a $MATRIX_LATTICES -rdisp 0.15 -sg $SG $NOSTORE $DR
