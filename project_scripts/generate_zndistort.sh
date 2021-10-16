../aiida_distort_structures.py -vs -0.01,-0.005,0.005,0.01 -in ZnEquilibrium_structures -out ZnEquilibrium-Distort_structures
../aiida_distort_structures.py -vs -0.01,-0.005,0.005,0.01 -in TphaseEquilib_structures -sg ZnEquilibrium-Distort_structures  -rdisp 0.15 -nrs 5 --max_atoms 150
