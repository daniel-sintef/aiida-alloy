LATTICE=4.04
MATRIX=Al
../aiida_create_surface_structures.py -a $LATTICE -me $MATRIX -l_surf FCC_111 -pzr 3 -vact 8 -sg Surface_structures
../aiida_create_surface_structures.py -a $LATTICE -me $MATRIX -l_surf FCC_110 -pzr 5 -vact 8 -sg Surface_structures
../aiida_create_surface_structures.py -a $LATTICE -me $MATRIX -l_surf FCC_100 -pzr 3 -vact 8 -sg Surface_structures

#LATTICE=3.62
#MATRIX=Cu
#../aiida_create_surface_structures.py -a $LATTICE -me $MATRIX -l_surf FCC_111 -pzr 3 -vact 0,6,0.25 -sg Surface_structures
#../aiida_create_surface_structures.py -a $LATTICE -me $MATRIX -l_surf FCC_110 -pzr 5 -vact 0,6,0.25 -sg Surface_structures
#../aiida_create_surface_structures.py -a $LATTICE -me $MATRIX -l_surf FCC_100 -pzr 3 -vact 0,6,0.25 -sg Surface_structures

#LATTICE=3.62
#MATRIX=Cu
#for PZR in {3..6}; do
#  ../aiida_create_surface_structures.py -a $LATTICE -me $MATRIX -l_surf FCC_111 -pzr $PZR -vact 8 -sg Surface_structures
#  ../aiida_create_surface_structures.py -a $LATTICE -me $MATRIX -l_surf FCC_110 -pzr $PZR -vact 8 -sg Surface_structures
#  ../aiida_create_surface_structures.py -a $LATTICE -me $MATRIX -l_surf FCC_100 -pzr $PZR -vact 8 -sg Surface_structures
#  done

