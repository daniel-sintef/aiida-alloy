#!/usr/bin/env python
import aiida
aiida.try_load_dbenv()

import click
import ase
import glob
import json
import os
from aiida_create_solutesupercell_structures import *

@click.command()
@click.option('-pbd', '--phonopy_basedir', required=True,
               help="path to a directory containing a dump of OQMD entries")
@click.option('-sg', '--structure_group_name', required=True,
              help="Output AiiDA group to store created structures")
@click.option('-sgd', '--structure_group_description', default="",
              help="Description for output AiiDA group")
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints structures and extras but does not store anything")
def launch(phonopy_basedir, structure_group_name, structure_group_description, dryrun):
    """
    Load an 'OQMD' dump (created by an custom script). Expects a directory containing a set
    of OQMD_<ID> vasp-formatted POSCAR, and for each of these a corresponding OQMD_<ID>.json
    json file wich contains a dump of the meta-data
    """

    print "loading dataset: {} to group: {}".format(phonopy_basedir, structure_group_name)
    # Setup/Retrieve the Group
    structure_group = Group.get_or_create(name=structure_group_name,
                            description=structure_group_description)[0]

    #Loop all relevant phonopy dirs
    phonopydirs = glob.glob(phonopy_basedir+'/PHONOPY_*')
    for phonopydir in phonopydirs:
        phonopy_json = os.path.join(phonopydir, "AiiDA.json")
        if not os.path.isfile(phonopy_json):
            print("No meta file found for {} refusing to load".format(
                   phonopy_dumpfile))
            continue
        with open(phonopy_json, 'r') as fp:
            phonopy_meta = json.load(fp)
        #Loop all relevant phonopy displacement poscar
        phonopy_disp_poscars = glob.glob(phonopydir+'/POSCAR-*')
        for phonopy_disp_poscar in phonopy_disp_poscars:
            phonopy_structure = ase.io.read(phonopy_disp_poscar, format='vasp')
            phonopy_disp_id = int(os.path.basename(phonopy_disp_poscar).split('-')[-1])
            #storing the phonopy displacement id, corresponding with phonopy_disp.yaml
            phonopy_meta['phonopy_disp_id'] = phonopy_disp_id
            phonopy_meta['phonopy_dirname'] = os.path.dirname(phonopy_disp_poscar)
            store_asestructure(phonopy_structure, phonopy_meta, structure_group, dryrun)

if __name__ == "__main__":
    launch()
