#!/usr/bin/env python
import aiida
aiida.try_load_dbenv()
from aiida.orm import Group
from aiida_create_solutesupercell_structures import *
import ase
import ase.build
import click
import numpy as np

def get_displacements_array(displacement):
    if len(displacement.split(',')) == 3:
       d_min,d_max,d_inc = displacement.split(',')
       displacements = np.arange(
                           float(d_min),float(d_max),float(d_inc)
                           )
    else:
       displacements = np.array([float(displacement)])
    return displacements


@click.command()
@click.option('-bsz', '--box_size', default=15,
              help="Length of (cubic) box in Ang.")
@click.option('-dsep', '--dimer_separation', required=True,
              help="The separation between the dimers in Ang"
                   "The notation is: start,end,increment or displacement_value.")
@click.option('-fse', '--firstdimer_elements', required=True,
              help="The first dimer."
              " Can pass a list of elements using comma seperation"
              " E.g. Mg,Si,Cu.")
@click.option('-sse', '--seconddimer_elements', required=True,
              help="Second dimer element, created at unique distances from the first dimer."
              " Can pass a list of elements using comma seperation"
              " E.g. Mg,Si,Cu. Can specify the creation of a vacancy using 'Vac'"
              " NOTE: will not generate symmetrically equivalent structures."
              " E.g. if Mg-Si dimer has been generated the script will skip Si-Mg")
@click.option('-sg', '--structure_group_name', required=True,
              help="Output AiiDA group to store created structures")
@click.option('-sgd', '--structure_group_description', default="",
              help="Description for output AiiDA group")
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints structures and extras but does not store anything")
def launch(box_size, dimer_separation,
           firstdimer_elements, seconddimer_elements,
           structure_group_name, structure_group_description,
           dryrun):
    """
    Script for creating surface structures for a given size and matrix element. Generates
    a set of structures with varying vacuum thickness
    """
    if not dryrun:
        structure_group = Group.get_or_create(
                             name=structure_group_name, description=structure_group_description)[0]
    else:
        structure_group = None

    box_size = float(box_size)
    extras = {'box_size':box_size}
    base_structure = ase.Atoms(cell=([box_size, box_size, box_size]), pbc=True)

    seperation_distances = get_displacements_array(dimer_separation)

    firstdimer_elements = prep_elementlist(firstdimer_elements)
    seconddimer_elements = prep_elementlist(seconddimer_elements)

    previous_firstdimer_elements = [] # to avoid duplication in generated structures
    for firstdimer_element in firstdimer_elements:
        dimer_structure = base_structure.copy()

        # create and store single atom version
        first_atom = ase.Atom(firstdimer_element)
        dimer_structure.append(first_atom)
        extras['firstdimer_element'] = firstdimer_element
        extras['second_element'] =  ""
        extras['atomatom_distance'] = ""
        store_asestructure(dimer_structure, extras,
                           structure_group, dryrun)

        dimer_structure.append(first_atom) #This will be changed later
        for seconddimer_element in seconddimer_elements:
            dimer_structure[1].symbol = seconddimer_element

            if seconddimer_element in previous_firstdimer_elements:
                continue
            extras['second_element'] =  seconddimer_element
            for distance in seperation_distances:
                dimer_structure[1].position[0] = distance
                extras['atomatom_distance'] = distance
                store_asestructure(dimer_structure, extras,
                                   structure_group, dryrun)

        previous_firstdimer_elements += [firstdimer_element]

if __name__ == "__main__":
   launch()
