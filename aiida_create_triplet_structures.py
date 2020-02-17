#!/usr/bin/env python
import aiida
aiida.try_load_dbenv()
from aiida.orm import Group
from aiida_create_solutesupercell_structures import *
import ase
import ase.build
import click
import numpy as np
import random
import itertools


@click.command()
@click.option('-a', '--lattice_size', required=True,
              help="lattice length (in Ang) to use")
@click.option('-spsh', '--supercell_shape', required=True,
              help="shape of the supercell to use, format: Nx,Ny,Nz")
@click.option('-me', '--matrix_element', required=True,
              help="element to be used as the matrix")
@click.option('-te', '--triplet_elements', required=True,
              help="elements to be used for the solute triplets")
@click.option('-sg', '--structure_group_name', required=True,
              help="Output AiiDA group to store created structures")
@click.option('-sgd', '--structure_group_description', default="",
              help="Description for output AiiDA group")
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints structures and extras but does not store anything")
def launch(lattice_size, supercell_shape, matrix_element,
           triplet_elements,
           structure_group_name, structure_group_description,
           dryrun):
    """
    Script for generating solute triplets 
    """
    if not dryrun:
        structure_group = Group.get_or_create(
                             name=structure_group_name, description=structure_group_description)[0]
    else:
        structure_group = None

    lattice_size = float(lattice_size)

    supercell_shape = supercell_shape.split(',')
    if len(supercell_shape) != 3:
        sys.exit("supercell_shape must be of the form Nx,Ny,Nz")

    base_extras = {
        'lattice_size':lattice_size,
        'supercell_shape':supercell_shape,
        'matrix_element':matrix_element,
                  }

    pure_structure = gen_ase_supercell(lattice_size, supercell_shape, matrix_element)
    pure_extras = copy.deepcopy(base_extras)

    triplet_elements = prep_elementlist(triplet_elements)
    triplets = list(itertools.combinations_with_replacement(
                      triplet_elements, 3))
    for triplet in triplets:
        triplet_extras = copy.deepcopy(base_extras)
        triplet_extras['triplet'] = triplet

        triplet_structure = copy.deepcopy(pure_structure)
        triplet_structure[0].symbol = triplet[0]
        triplet_structure[1].symbol = triplet[1]
        triplet_structure[2].symbol = triplet[2]

        store_asestructure(triplet_structure, triplet_extras,
                           structure_group, dryrun)



if __name__ == "__main__":
   launch()
