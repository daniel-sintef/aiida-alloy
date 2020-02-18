#!/usr/bin/env python
import click
import ase
import ase.build
import json
import os
import aiida
aiida.try_load_dbenv()
from aiida.orm import QueryBuilder
from aiida.orm import Node, Group
from aiida.orm import load_node

def get_structurenode_metadict(structure_node):
    meta_dict = structure_node.get_extras()
    meta_dict['uuid'] = structure_node.uuid
    return meta_dict


def export_structure_node(structure_node, output_path, ase_format='vasp'):
    ase_structure = structure_node.get_ase()
    meta_dict = get_structurenode_metadict(structure_node)

    ase_structure.write(output_path, format=ase_format)

    with open(output_path+'.json', 'w') as fp:
        json.dump(meta_dict, fp)
    return

def get_allnodes_fromgroup(group_name):
    qb = QueryBuilder()
    qb.append(Group, filters={'name': group_name}, tag='g')
    qb.append(Node, tag='job', with_group='g')
    all_nodes = [x[0] for x in qb.all()]
    return all_nodes

@click.command()
@click.option('-od', '--output_dir', required=True
         , help="Directory to dump output files")
@click.option('-gn', '--group_name',default=None,
              type=str, help="Group to export identified by name")
@click.option('-u', '--uuid', default=None,
              type=str, help="Structure to export by uuid")
def createjob(output_dir, group_name, uuid):
    '''
    Dumps the contents of an AiiDA group into a directory. Creates a set of
    VASP-poscar files with the name AIIDA_<AiiDA-pk> for each of
    these coordinate files, a corresponding .json file is created
    which contains the meta data for the entry.
    '''
    try:
        os.makedirs(output_dir)
    except Exception:
        pass

    if group_name is not None:
        all_entries = get_allnodes_fromgroup(group_name)
    elif uuid:
        all_entries = [load_node(uuid)]
    else:
        raise Exception("You must provide either group_name or uuid")

    for structure_node in all_entries:
        base_filename = "AIIDA_{}".format(structure_node.pk)
        output_path = os.path.join(output_dir, base_filename)
        if os.path.isfile(output_path):
           print(("{} already exists skipping".format(output_path)))
           continue

        try:
            export_structure_node(structure_node, output_path)
        except Exception:
            print(("Failed to parse {}".format(structure_node)))
            fail_dict = get_metadict_from_oqmdentry(structure_node)
            with open(output_path+'.EXPORT_FAIL', 'w') as fp:
                json.dump(fail_dict, fp)


    return


if __name__ == "__main__":
    createjob()
