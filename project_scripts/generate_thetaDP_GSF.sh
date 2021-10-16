STRUCTURE_UUID='2dcf7617-cd0f-4f19-a54e-5533c5af3798'
GROUPNAME=ThetaDoublePrime_GSF_Structures
GENSTRUCTURE_CMD="../aiida_create_stackingfault_structures.py -sg $GROUPNAME -cstr  $STRUCTURE_UUID"

#DRYRUN=--dryrun
DRYRUN=
PERIODIC_Z="--periodic_zrepeats 8"
for i in {0..3}; do
  DX=$(echo "scale=8; 0" | bc)
  DY=$(echo "scale=8; 0+$i/4" | bc)
  $GENSTRUCTURE_CMD  -dx $DX -dy $DY $DRYRUN $PERIODIC_Z

  DX=$(echo "scale=8; 1/3" | bc)
  DY=$(echo "scale=8; 1/12+$i/4" | bc)
  $GENSTRUCTURE_CMD  -dx $DX -dy $DY $DRYRUN $PERIODIC_Z

  DX=$(echo "scale=8; 2/3" | bc)
  DY=$(echo "scale=8; 2/12+$i/4" | bc)
  $GENSTRUCTURE_CMD  -dx $DX -dy $DY $DRYRUN $PERIODIC_Z
done

DRYRUN=
PERIODIC_Z="--periodic_zrepeats 12"
for i in {0..3}; do
  DX=$(echo "scale=8; 0" | bc)
  DY=$(echo "scale=8; 0+$i/4" | bc)
  $GENSTRUCTURE_CMD  -dx $DX -dy $DY $DRYRUN $PERIODIC_Z

  DX=$(echo "scale=8; 1/3" | bc)
  DY=$(echo "scale=8; 1/12+$i/4" | bc)
  $GENSTRUCTURE_CMD  -dx $DX -dy $DY $DRYRUN $PERIODIC_Z

  DX=$(echo "scale=8; 2/3" | bc)
  DY=$(echo "scale=8; 2/12+$i/4" | bc)
  $GENSTRUCTURE_CMD  -dx $DX -dy $DY $DRYRUN $PERIODIC_Z
done
