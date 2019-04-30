#!/usr/bin/env python
import aiida
aiida.try_load_dbenv()

import click
import ase
import glob
import json
import os
from aiida_create_solutesupercell_structures import *

def add_parentstructure_extras(structurenode, parent_uuid):
    # NOTE: consider adding a check if parent_extras is already assigned
    structure_extras = structurenode.get_extras()
    parent_extras = load_node(parent_uuid).get_extras()
    for key, value in parent_extras.items():
        if key not in structure_extras:
            structurenode.set_extra(key, value)
    structurenode.set_extra('parent_extras', True)
    return


@click.command()
@click.option('-od', '--oqmd_dumpdir', required=True,
               help="path to a directory containing a dump of OQMD entries")
@click.option('-sg', '--structure_group_name', required=True,
              help="Output AiiDA group to store created structures")
@click.option('-sgd', '--structure_group_description', default="",
              help="Description for output AiiDA group")
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints structures and extras but does not store anything")
def launch(oqmd_dumpdir, structure_group_name, structure_group_description, dryrun):
    """
    Load an 'OQMD' dump (created by an custom script). Expects a directory containing a set
    of OQMD_<ID> vasp-formatted POSCAR, and for each of these a corresponding OQMD_<ID>.json
    json file wich contains a dump of the meta-data
    """

    print "loading dataset: {} to group: {}".format(oqmd_dumpdir, structure_group_name)
    # Setup/Retrieve the Group
    structure_group = Group.get_or_create(name=structure_group_name,
                            description=structure_group_description)[0]

    oqmd_dumpfiles = glob.glob(oqmd_dumpdir+'/*')
    oqmd_dumpfiles = filter(lambda x: '.' not in os.path.basename(x), oqmd_dumpfiles)
    print oqmd_dumpfiles
    for oqmd_dumpfile in oqmd_dumpfiles:
        oqmd_structure = ase.io.read(oqmd_dumpfile, format='vasp')

        oqmd_json = oqmd_dumpfile+'.json'
        if not os.path.isfile(oqmd_json):
            print("No meta file found for {} refusing to load".format(
                   oqmd_dumpfile))
            continue
        with open(oqmd_json, 'r') as fp:
            oqmd_meta = json.load(fp)

        store_asestructure(oqmd_structure, oqmd_meta, structure_group, dryrun)

if __name__ == "__main__":
    launch()
