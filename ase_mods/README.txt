in order for load_runner_dataset.py to work, one needs to modify the ase install to be able
to read 'runner' data

1) copy the runner.py format to the io dir of ase
$cp ase_mods/runner.py ${VIRTUAL_ENV}/lib/python2.7/site-packages/ase/io/

2) add the 'runner_line' to the list of formats in the all_fromats file in
   $VIRTUAL_ENV/lib/python2.7/site-packages/ase/io/formats.py

runner_line: 'runner': ('Runner input file', '+F'),

ase_mods/runner.py is provided as an example, DO NOT copy this file, merely
modify the existing version to include the 'runner_line'

