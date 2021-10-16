#!/usr/bin/env python
import aiida
aiida.load_profile()

import click
import ase
import glob
import json
import os
from aiida_create_solutesupercell_structures import *

@click.command()
@click.option('-od', '--oqmd_dumpdir', required=True,
               help="path to a directory containing a dump of OQMD entries")
@click.option('-e', '--extras', required=False,
              help="Add extras, each key,label is joined by a comman and seperated"
                   " by pipes. e.g. key1,label1|key2,label2")
@click.option('-sg', '--structure_group_label', required=True,
              help="Output AiiDA group to store created structures")
@click.option('-sgd', '--structure_group_description', default="",
              help="Description for output AiiDA group")
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints structures and extras but does not store anything")
def launch(oqmd_dumpdir, structure_group_label, structure_group_description, extras,  dryrun):
    """
    Load an 'OQMD' dump (created by an custom script). Expects a directory containing a set
    of OQMD_<ID> vasp-formatted POSCAR, and for each of these a corresponding OQMD_<ID>.json
    json file wich contains a dump of the meta-data
    """
    if extras is not None:
        extras = {y[0]:y[1] for y in [x.split(',') for x in extras.split('|')]}
    else:
        extras = {}

    print("loading dataset: {} to group: {}".format(oqmd_dumpdir, structure_group_label))
    # Setup/Retrieve the Group
    structure_group = Group.objects.get_or_create(label=structure_group_label,
                            description=structure_group_description)[0]

    oqmd_dumpfiles = glob.glob(oqmd_dumpdir+'/*')
    oqmd_dumpfiles = [x for x in oqmd_dumpfiles if '.' not in os.path.basename(x)]
    print(oqmd_dumpfiles)
    for oqmd_dumpfile in oqmd_dumpfiles:
        oqmd_structure = ase.io.read(oqmd_dumpfile, format='vasp')

        oqmd_json = oqmd_dumpfile+'.json'
        if not os.path.isfile(oqmd_json):
            print(("No meta file found for {} refusing to load".format(
                   oqmd_dumpfile)))
            continue
        with open(oqmd_json, 'r') as fp:
            oqmd_meta = json.load(fp)

        for k,v in extras.items():
            if k not in oqmd_meta:
                oqmd_meta[k] = v
        store_asestructure(oqmd_structure, oqmd_meta, structure_group, dryrun)

if __name__ == "__main__":
    launch()
