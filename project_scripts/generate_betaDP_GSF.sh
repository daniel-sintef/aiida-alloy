GROUPNAME=BetaDPGSF_structures
GENSTRUCTURE_CMD="../aiida_create_stackingfault_structures.py -sg $GROUPNAME"

#DRYRUN=
PERIODIC_Z="--periodic_zrepeats 8"
STRUCTURE_UUID=a3148955-c626-4c66-984e-cd4ea056f316
STRUCTURE_CMD="-cstn $STRUCTURE_UUID"
DISPLACEMENTS="0.0,0.0
0.25,0.00
0.50,0.00
0.00,0.40
0.30,0.40
0.50,0.60
0.15,0.40
0.60,0.40
0.20,0.60
"
for DISPLACEMENT in $DISPLACEMENTS; do
  DX=$(echo $DISPLACEMENT | awk -F, '{print $1}')
  DY=$(echo $DISPLACEMENT | awk -F, '{print $2}')
  $GENSTRUCTURE_CMD  $DRYRUN -dx $DX -dy $DY $PERIODIC_Z $STRUCTURE_CMD
done
