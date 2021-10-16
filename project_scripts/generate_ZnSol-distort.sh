#DRYRUN="--dryrun"
DRYRUN=""

NORMSTRAIN=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.010,0.010,6), decimals=5)]))")
SHEARSTRAIN=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(0.000,0.001,6), decimals=5)]))")
../aiida_distort_structures.py --input_structures  '5e2faa1f-7e0b-465b-91cc-2dd6f5bd8b03' \
  $DRYRUN \
  --norm_strains $NORMSTRAIN \
  --shear_strains $SHEARSTRAIN \
  -sg DiluteZnInAl-distorted_structures \
  --max_atoms 109 \
  $DRYRUN

VOLRANGE=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.005,0.01,3), decimals=5)]))")
../aiida_distort_structures.py --input_structures  '5e2faa1f-7e0b-465b-91cc-2dd6f5bd8b03' \
   --volumetric_strains $VOLRANGE \
  -sg DiluteZnInAl-distorted_structures \
  --max_atoms 109 \
  $DRYRUN
../aiida_distort_structures.py --input_structures  '5e2faa1f-7e0b-465b-91cc-2dd6f5bd8b03' \
   --volumetric_strains $VOLRANGE \
  -rdisp 0.05 -nrs 5 \
  -sg DiluteZnInAl-distorted_structures \
  --max_atoms 109 \
  $DRYRUN

