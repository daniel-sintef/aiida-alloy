# Testing Al stacking fault energy

#-sg SoluteStackingFault_structures
../aiida_create_stackingfault_structures.py \
  -a 4.04 -l_surf FCC_111 -me Al -se Zn \
  -pxr 4 -pyr 4 -pzr 3 \
  --primitive \
  -sg SoluteStackingFault_structures

../aiida_create_stackingfault_structures.py \
  -a 4.04 -l_surf FCC_111 -me Al -se Zn \
  -pxr 4 -pyr 4 -pzr 3 \
  --primitive --refsolute \
  -sg SoluteStackingFault_structures

##-sg SoluteStackingFault_structures
#../aiida_create_stackingfault_structures.py \
#  -a 4.04 -l_surf FCC_111 -me Al -se Cu,Mg,Si \
#  -pxr 4 -pyr 4 -pzr 3 \
#  --primitive \
#  -sg SoluteStackingFault_structures
#
#../aiida_create_stackingfault_structures.py \
#  -a 4.04 -l_surf FCC_111 -me Al -se Cu,Mg,Si \
#  -pxr 4 -pyr 4 -pzr 3 \
#  --primitive --refsolute \
#  -sg SoluteStackingFault_structures

#../aiida_create_stackingfault_structures.py \
#  -a 3.62 -l_surf FCC_111 -me Cu -se Al \
#  -pxr 3 -pyr 3 -pzr 3 \
#  --primitive --testsolute_layer \
#  -sg SoluteStackingFault_structuresCu
#
#../aiida_create_stackingfault_structures.py \
#  -a 3.62 -l_surf FCC_111 -me Cu -se Al \
#  -pxr 3 -pyr 3 -pzr 3 \
#  --primitive  \
#  -sg SoluteStackingFault_structuresCu

##out-plane benchmarking
#for Z_REPEAT in {2..5}; do
#  ../aiida_create_stackingfault_structures.py \
#  -a 4.04 -l_surf FCC_111 -me Al -se Cu,Mg,Si,Vac \
#  -pxr 4 -pyr 4 -pzr $Z_REPEAT --primitive \
#  -sg SoluteStackingFault_structures \
#  --testsolute_layer
#
#  done


## TESTING in-plane thickness
#MAX_SOLUTE_LAYER=0
#for XY_REPEATS in {2..4}; do
#
#  ../aiida_create_stackingfault_structures.py \
#  -a 4.04 -l_surf FCC_111 -me Al -se Cu,Mg,Si,Vac \
#  -pxr $XY_REPEATS -pyr $XY_REPEATS -pzr 2 \
#  --primitive --maxsolute_layer $MAX_SOLUTE_LAYER \
#  -sg SoluteStackingFault_structures
#
#  ../aiida_create_stackingfault_structures.py \
#  -a 4.04 -l_surf FCC_111 -me Al -se Cu,Mg,Si,Vac \
#  -pxr $XY_REPEATS -pyr $XY_REPEATS -pzr 2 \
#  --primitive --refsolute \
#  -sg SoluteStackingFault_structures
#
#  ../aiida_create_stackingfault_structures.py \
#  -a 3.62 -l_surf FCC_111 -me Cu -se Al \
#  -pxr $XY_REPEATS -pyr $XY_REPEATS -pzr 2 \
#  --primitive --maxsolute_layer $MAX_SOLUTE_LAYER \
#  -sg SoluteStackingFault_structures
#
#  ../aiida_create_stackingfault_structures.py \
#  -a 3.62 -l_surf FCC_111 -me Cu -se Al \
#  -pxr $XY_REPEATS -pyr $XY_REPEATS -pzr 2 \
#  --primitive --refsolute \
#  -sg SoluteStackingFault_structures
#  done
#
