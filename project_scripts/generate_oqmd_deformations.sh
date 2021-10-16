#Structure uuid
#b6707b66-933b-42d9-af06-0d411742b028 Mg
#2854624b-7ddf-4968-b02f-7b2e54777859 Si
#593927fd-6816-46ac-b585-ba5051b6bdaf Al
#b97f4359-860b-411b-adb3-081e357fbe9c Cu
#INPUT_STRUCTURES="b6707b66-933b-42d9-af06-0d411742b028,2854624b-7ddf-4968-b02f-7b2e54777859,593927fd-6816-46ac-b585-ba5051b6bdaf,b97f4359-860b-411b-adb3-081e357fbe9c,e538fb62-69d2-49b6-9c9d-9503938298c5,e708f2d1-ac01-4f8a-813c-c3639caf631b,a38c8be7-e52f-4e13-92c6-7866a4c6b294,09542054-df14-4799-a1ab-f1e424b154c8,63d867de-499a-4265-9dae-b8099f43adac"
# Just the pure elements
INPUT_STRUCTURES="b6707b66-933b-42d9-af06-0d411742b028,2854624b-7ddf-4968-b02f-7b2e54777859,593927fd-6816-46ac-b585-ba5051b6bdaf,b97f4359-860b-411b-adb3-081e357fbe9c"

#DRYRUN="--dryrun"
DRYRUN=""

# Lot's of Al, finely-tuned around the elastic constants
NORMSTRAIN=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.01,0.01,50), decimals=5)]))")
SHEARSTRAIN=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(0.000,0.01,50), decimals=5)]))")
echo $SHEARSTRAIN
../aiida_distort_structures.py --input_structures 593927fd-6816-46ac-b585-ba5051b6bdaf \
  $DRYRUN \
  -sg EOS-SmallStrainOnly_structures \
  --norm_strains $NORMSTRAIN \
  --shear_strains $SHEARSTRAIN \
  --use_conventional_structure \

## Lot's of Al, finely-tuned around the elastic constants
#VOLRANGE=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.005,0.005,10), decimals=5)]))")
#NORMSTRAIN=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.01,0.01,40), decimals=5)]))")
#SHEARSTRAIN=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(0.000,0.01,40), decimals=5)]))")
#../aiida_distort_structures.py --input_structures 593927fd-6816-46ac-b585-ba5051b6bdaf \
#  $DRYRUN \
#  -sg EOS-SmallStrainOnly_structures \
#  --norm_strains $NORMSTRAIN \
#  --shear_strains $SHEARSTRAIN \
#  --use_conventional_structure \
#  --volumetric_strains $VOLRANGE \

#VOLRANGE=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.5,0.5,20), decimals=5)]))")
#../aiida_distort_structures.py  \
#   --volumetric_strains $VOLRANGE \
#  --repeat_expansion 2,2,2 \
#  --input_structures $INPUT_STRUCTURES \
#  -sg OQMD-Distorted_structures \
#  $DRYRUN
#
##Pure elastic
#DEFORMRANGE=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.4,0.4,12), decimals=5)]))")
#../aiida_distort_structures.py --input_structures $INPUT_STRUCTURES \
#  --repeat_expansion 2,2,2 \
#  -sg OQMD-Distorted_structures  --structure_comments "pure elastic" \
#  --elastic_strains $DEFORMRANGE \
#  $DRYRUN

##Combined elastic, small amount
#VOLRANGE=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.005,0.005,3), decimals=5)]))")
#DEFORMRANGE=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.005,0.005,3), decimals=5)]))")
#./aiida_distort_structures.py --input_structures $INPUT_STRUCTURES \
#  -sg OQMD-Distorted_structures -rdisp 0.05 -nrs 4 \
#  --elastic_strains $DEFORMRANGE --volumetric_strains $VOLRANGE \
#  --structure_comments "combined everything" \
#  $DRYRUN

### Pure (extreme) volume
#VOLRANGE=$(python -c "import numpy as np;print(','.join([str(x) for x in np.round(np.linspace(-0.05,0.1,20), decimals=5)]))")
#./aiida_distort_structures.py --input_structures $INPUT_STRUCTURES \
#    -sg OQMD-Distorted_structures \
#    --volumetric_strains $VOLRANGE \
#    --structure_comments "pure volume" \
#    $DRYRUN
#
##Large(r) MD simulation
#./aiida_distort_structures.py --input_structures $INPUT_STRUCTURES \
#    -sg OQMD-Distorted_structures -rdisp 0.05 -nrs 5 --repeat_expansion 2,2,2 \
#    --structure_comments "random-random" \
#    $DRYRUN
#./aiida_distort_structures.py --input_structures $INPUT_STRUCTURES \
#    -sg OQMD-Distorted_structures -rdisp 0.10 -nrs 5 --repeat_expansion 2,2,2 \
#    --structure_comments "random-random" \
#    $DRYRUN
#./aiida_distort_structures.py --input_structures $INPUT_STRUCTURES \
#    -sg OQMD-Distorted_structures -rdisp 0.15 -nrs 5 --repeat_expansion 2,2,2 \
#    --structure_comments "random-random" \
#    $DRYRUN
#./aiida_distort_structures.py --input_structures $INPUT_STRUCTURES \
#    -sg OQMD-Distorted_structures -rdisp 0.20 -nrs 5 --repeat_expansion 2,2,2 \
#    --structure_comments "random-random" \
#    $DRYRUN
