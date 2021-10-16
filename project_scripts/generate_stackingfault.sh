# Testing Al stacking fault energy
#../aiida_create_stackingfault_structures.py --lattice_size 4.04 --matrix_element Al --lattice_and_surface FCC_111 --periodic_zrepeats 5 --special_pointsonly --structure_group_name StackingFault_structures
#../aiida_create_stackingfault_structures.py -a 3.62 -me Cu -l_surf FCC_111 -pzr 5 -dy 0,1,0.01 --structure_group_name StackingFault_structures
../aiida_create_stackingfault_structures.py -a 4.04 -me Al -l_surf FCC_111 -pzr 5 --special_pointsonly --structure_group_name StackingFault_structures
