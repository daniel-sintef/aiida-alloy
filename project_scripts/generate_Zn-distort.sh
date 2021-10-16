#DRYRUN="--dryrun"
DRYRUN=""

NORMSTRAIN=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.025,0.025,15), decimals=5)]))")
SHEARSTRAIN=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(0.000,0.025,15), decimals=5)]))")
echo $SHEARSTRAIN
../aiida_distort_structures.py --input_structures 218c4f5c-8b6c-4981-96f2-ea9b9e1940d5 \
  $DRYRUN \
  -sg Zn-distortion_structures \
  --norm_strains $NORMSTRAIN \
  --shear_strains $SHEARSTRAIN \
  --use_conventional_structure \
  $DRYRUN

VOLRANGE=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.02,0.02,10), decimals=5)]))")
../aiida_distort_structures.py  \
   --volumetric_strains $VOLRANGE \
  --repeat_expansion 2,2,2 \
  --input_structures 218c4f5c-8b6c-4981-96f2-ea9b9e1940d5 \
  -rdisp 0.05 -nrs 5 --repeat_expansion 3,2,2 \
  -sg Zn-distortion_structures \
  $DRYRUN

