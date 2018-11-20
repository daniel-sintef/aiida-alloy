#!/usr/bin/env python
import aiida
aiida.load_dbenv()
from aiida.orm.group import Group
from aiida.orm.data.structure import StructureData
import ase
import ase.io
import click


@click.command()
@click.option('-d', '--dataset_path', required=True)
@click.option('-gn', '--group_name', required=True)
@click.option('-gd', '--group_description', default="")
def launch(dataset_path, group_name, group_description):
    print "loading dataset: {} to group: {}".format(dataset_path, group_name)

    # Setup/Retrieve the Group
    g = Group.get_or_create(name=group_name, description=group_description)[0]
    # Loop over structures in the dataset_path, storing the nodes then adding them to the group
    i = 0
    while True:
        try:
            ase_structure = ase.io.read(dataset_path, index=i, format="runner")
        except StopIteration:
            break
        # setup and store the ase structure as an aiida StructureData node
        aiida_structure = StructureData()
        aiida_structure.set_ase(ase_structure)
        aiida_structure_stored = aiida_structure.store()

        # add in the dataset_path line if possible
        try:
            structure_path = ase_structure.comment.strip().split()[-1][3:]
            aiida_structure_stored.set_extra("structure_path", structure_path)
        except AttributeError:
            print "could not set structure_path on {}".format(ase_structure)
            pass

        # add in the chemical formula and number of atoms if possible
        try:
            aiida_structure_stored.set_extra("num_atoms",
                                             len(ase_structure))
            aiida_structure_stored.set_extra("chem_formula",
                          ase_structure.get_chemical_formula())
        except AttributeError:
            print "could not set either num_atoms or chemical_formula " \
                  " on {}".format(ase_structure)
            pass

        # add the structure to the group
        g.add_nodes(aiida_structure_stored)
        g.store()
        i += 1


if __name__ == "__main__":
    try:
        ase.io.read("", format="runner")
    except ValueError:
        raise ValueError("You need a version of ase that can read runner files")
    except IOError:
        pass
    launch()
