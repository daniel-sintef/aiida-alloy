#!/usr/bin/env python
import aiida
aiida.load_dbenv()
from aiida.orm.group import Group
from aiida.orm.data.structure import StructureData
import ase
from ase.build import sort
import click
import copy
import numpy as np
import pandas as pd
import sys

'''
ase requires that all assigned symbols be on the periodic table. Here we have chosen
the highly unstable Nobelium (element 102) to internally represent vacancies.
'''
VACANCY_INTERNAL_SYMBOL="No"
VACANCY_USER_SYMBOL="Vac"


def gen_ase_supercell(lattice_size, supercell_shape, matrix_element):
    a1 = np.array([lattice_size,0.,0.])
    a2 = np.array([0.,lattice_size,0.])
    a3 = np.array([0.,0.,lattice_size])
    t1 = np.array([0.0,0.0,0.0])
    t2 = np.array([0.0,0.5,0.5])
    t3 = np.array([0.5,0.0,0.5])
    t4 = np.array([0.5,0.5,0.0])

    atoms = ase.Atoms('{}4'.format(matrix_element),
                      cell=[a1,a2,a3], pbc=True,
                      scaled_positions=[t1,t2,t3,t4])
    supercell_atoms = atoms*np.array(supercell_shape, dtype=int)
    return supercell_atoms

def return_nn_distanceAndIndex(ase_supercell):
    supercell_frame = pd.DataFrame(ase_supercell.get_positions(),
                                   columns=['x','y','z'])
    supercell_frame['distance'] = supercell_frame.apply(
                                    np.linalg.norm, axis=1)
    supercell_frame['distance'] = supercell_frame['distance'].round(9)

    supercell_shape = ase_supercell.get_cell()
    lx = supercell_shape[0][0]
    ly = supercell_shape[1][1]
    lz = supercell_shape[2][2]

    # we assume the first solute is always at the origin
    supercell_frame = supercell_frame[
                 (supercell_frame['x'] <= lx/2.) &
                 (supercell_frame['y'] <= ly/2.) &
                 (supercell_frame['z'] <= lz/2.) ]

    supercell_frame.sort_values('distance', inplace=True)
    supercell_frame.drop_duplicates(subset='distance',
                                    inplace=True)



    distance_values = supercell_frame['distance'].values
    index_values = supercell_frame.index.values

    nn_distanceindex_frame = pd.DataFrame(data=zip(distance_values,index_values),
                                          columns=['distances','indexes'])

    return nn_distanceindex_frame

def prep_elementlist(elementlist):
    elementlist = list(elementlist.split(','))
    elementlist = map(lambda x:x if x.lower()!= VACANCY_USER_SYMBOL.lower()
                      else VACANCY_INTERNAL_SYMBOL, elementlist)
    return elementlist

def store_asestructure(ase_structure, extras, structure_group, dryrun):
    ase_structure = sort(ase_structure)

    # convert any instances of vacancy internal symbol use back to user symbol use
    for key in extras:
        if extras[key] == VACANCY_INTERNAL_SYMBOL: extras[key] = VACANCY_USER_SYMBOL

    # delete all the vacancy sites prior to storage
    del ase_structure[[x.index for x in ase_structure
                       if x.symbol==VACANCY_INTERNAL_SYMBOL or
                          x.symbol==VACANCY_USER_SYMBOL]]

    if dryrun:
        print "structure: {}".format(ase_structure)
        print "extras: {}".format(extras)
    else:
        aiida_structure = StructureData()
        aiida_structure.set_ase(ase_structure)
        aiida_structure_stored = aiida_structure.store()
        for key in extras:
            aiida_structure_stored.set_extra(key, extras[key])

        aiida_structure_stored.set_extra("num_atoms", len(ase_structure))
        aiida_structure_stored.set_extra("chem_formula", ase_structure.get_chemical_formula())

        structure_group.add_nodes(aiida_structure_stored)

    return

@click.command()
@click.option('-a', '--lattice_size', required=True,
              help="lattice length (in Ang) to use")
@click.option('-spsh', '--supercell_shape', required=True,
              help="shape of the supercell to use, format: Nx,Ny,Nz")
@click.option('-me', '--matrix_element', required=True,
              help="element to be used as the ")
@click.option('-fse', '--firstsolute_elements', required=True,
              help="First solute element, always centered at the origin of the supercell."
              " Can pass a list of elements using comma seperation"
              " E.g. Mg,Si,Cu. Can specify the creation of a vacancy using 'Vac'")
@click.option('-sse', '--secondsolute_elements', required=True,
              help="Second solute element, created at unique distances from the first solute."
              " Can pass a list of elements using comma seperation"
              " E.g. Mg,Si,Cu. Can specify the creation of a vacancy using 'Vac'"
              " NOTE: will not generate symmetrically equivalent structures."
              " E.g. if Mg-Si solute solutes have been generated the script will skip Si-Mg")
@click.option('-sg', '--structure_group_name', required=True,
              help="Output AiiDA group to store created structures")
@click.option('-sgd', '--structure_group_description', default="",
              help="Description for output AiiDA group")
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints structures and extras but does not store anything")
def launch(lattice_size,
           supercell_shape, matrix_element,
           firstsolute_elements, secondsolute_elements,
           structure_group_name, structure_group_description, dryrun):

    lattice_size = float(lattice_size)

    supercell_shape = supercell_shape.split(',')
    if len(supercell_shape) != 3:
        sys.exit("supercell_shape must be of the form Nx,Ny,Nz")

    firstsolute_elements = prep_elementlist(firstsolute_elements)
    secondsolute_elements = prep_elementlist(secondsolute_elements)

    structure_group = Group.get_or_create(
                        name=structure_group_name, description=structure_group_description)[0]

    base_extras = {
        'lattice_size':lattice_size,
        'supercell_shape':supercell_shape,
        'matrix_element':matrix_element,
                  }

    pure_structure = gen_ase_supercell(lattice_size, supercell_shape, matrix_element)

    pure_extras = copy.deepcopy(base_extras)
    store_asestructure(pure_structure, pure_extras, structure_group, dryrun)

    nn_distanceindex_frame = return_nn_distanceAndIndex(pure_structure)

    previously_generated_firstsol_elements = [] # to avoid duplication in generated structures
    for firstsolute_element in firstsolute_elements:

        singlesol_structure = copy.deepcopy(pure_structure)
        singlesol_structure[0].symbol = firstsolute_element

        singlesol_extras = copy.deepcopy(pure_extras)
        singlesol_extras['sol1_element'] = firstsolute_element
        singlesol_extras['sol1_index'] = 0
        store_asestructure(singlesol_structure, singlesol_extras, structure_group, dryrun)

        for secondsolute_element in secondsolute_elements:
            # skip symmetrically equivalent structures
            if secondsolute_element in previously_generated_firstsol_elements:
                continue

            for i in range(1, len(nn_distanceindex_frame)):
                secondsol_structure = copy.deepcopy(singlesol_structure)
                secondsol_index = nn_distanceindex_frame['indexes'][i]
                secondsol_structure[secondsol_index].symbol = secondsolute_element

                secondsol_extras = copy.deepcopy(singlesol_extras)
                secondsol_extras['sol2_element'] = secondsolute_element
                secondsol_extras['sol2_index'] = secondsol_index
                secondsol_distance = nn_distanceindex_frame['distances'][i]
                secondsol_extras['sol2_distance'] = secondsol_distance
                store_asestructure(secondsol_structure, secondsol_extras,
                                   structure_group, dryrun)

        previously_generated_firstsol_elements += firstsolute_element

if __name__ == "__main__":
   launch()
