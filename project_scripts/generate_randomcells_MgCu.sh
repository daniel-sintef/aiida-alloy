MATRIX_ELEMENTS=Mg,Cu,Vac
MATRIX_LATTICES=4.58,3.65,-1
NUM_REPEATS=30
SG=Random-MgCu_structures
DR=

#no vacancy
echo "No vacancy"
CONCENTRATIONS="1.0,0.0,0.0
0.9,0.1,0.0
0.8,0.2,0.0
0.7,0.3,0.0
0.6,0.4,0.0
0.5,0.5,0.0
0.4,0.6,0.0
0.3,0.7,0.0
0.2,0.8,0.0
0.1,0.9,0.0
0.0,1.0,0.0"
for CONC in $CONCENTRATIONS; do
  ../aiida_create_randomsupercell_structures.py -me $MATRIX_ELEMENTS -a $MATRIX_LATTICES -c $CONC -rdisp 0.15 -ns $NUM_REPEATS -spsh 2,2,2 -sg $SG $DR
  done

#5% vacancy
echo "5% vacancy"
CONCENTRATIONS="0.95,0.0,0.05
0.875,0.075,0.05
0.775,0.175,0.05
0.675,0.275,0.05
0.575,0.375,0.05
0.475,0.475,0.05
0.375,0.575,0.05
0.275,0.675,0.05
0.175,0.775,0.05
0.075,0.875,0.05
0.0,0.95,0.05"
for CONC in $CONCENTRATIONS; do
  ../aiida_create_randomsupercell_structures.py -me $MATRIX_ELEMENTS -a $MATRIX_LATTICES -c $CONC -rdisp 0.15 -ns $NUM_REPEATS -spsh 2,2,2 -sg $SG $DR
  done

#10% vacancy
echo "10% vacancy"
CONCENTRATIONS="0.9,0.0,0.1
0.85,0.05,0.1
0.75,0.15,0.1
0.65,0.25,0.1
0.55,0.35,0.1
0.45,0.45,0.1
0.35,0.55,0.1
0.25,0.65,0.1
0.15,0.75,0.1
0.1,0.8,0.1
0.0,0.9,0.1"
for CONC in $CONCENTRATIONS; do
  ../aiida_create_randomsupercell_structures.py -me $MATRIX_ELEMENTS -a $MATRIX_LATTICES -c $CONC -rdisp 0.15 -ns $NUM_REPEATS -spsh 2,2,2 -sg $SG $DR
  done
