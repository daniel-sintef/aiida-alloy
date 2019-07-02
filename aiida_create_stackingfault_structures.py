#!/usr/bin/env python
import aiida
aiida.try_load_dbenv()
from aiida.orm.group import Group
from aiida.orm import load_node
from aiida_create_solutesupercell_structures import *
import ase
import ase.build
from ase.geometry import get_layers
import click
import numpy as np
import pandas as pd

def get_displacements_array(displacement):
    if len(displacement.split(',')) == 3:
       d_min,d_max,d_inc = displacement.split(',')
       if d_max > 1:
           print "WARNING: max displacement {} is larger than 1".format(displacement)
       displacements = np.arange(float(d_min),float(d_max),float(d_inc))
    else:
       displacements = np.array([float(displacement)])
    return displacements

def get_layer_frame(structure, miller_index):
    layer_indexes, layer_distances = get_layers(structure, miller_index)
    repeating_layer_distances = []
    for i in range(len(layer_indexes)):
        layer_distance = layer_distances[layer_indexes[i]]
        repeating_layer_distances.append(layer_distance)
    repeating_layer_distances = np.array(repeating_layer_distances)

    layer_dict = {'layer_index': layer_indexes,
                  'layer_distance': repeating_layer_distances}
    layer_frame = pd.DataFrame(data=layer_dict)
    layer_frame.index.name = "structure_index"
    return layer_frame


@click.command()
@click.option('-a', '--lattice_size', required=False,
              help="lattice length (in Ang) to use")
@click.option('-me', '--matrix_element', required=False,
              help="element to be used as the matrix")
@click.option('-l_surf', '--lattice_and_surface',
              type=click.Choice(["FCC_111"]),
              help="lattice and surface to be used. "
              "FCC_111: <112>(x) <110>(y) <111>(z)", required=False)
@click.option('-cstr', '--custom_structure',
              help="Node containing a custom structure. "
              "Must have x_direction, y_direction and surface_plane specified in the extras",
              required=False)
@click.option('-pxr', '--periodic_xrepeats', default=1, type=int,
              help="periodic repeats in the x (a1) direction")
@click.option('-pyr', '--periodic_yrepeats', default=1, type=int,
              help="periodic repeats in the y (a2) direction")
@click.option('-pzr', '--periodic_zrepeats', default=1, type=int,
              help="periodic repeats in the z (a3) direction")
@click.option('-dx', '--displacement_x', default='0',
              help="Displacements to be made in the x (a1) direction. "
                   "The amount is expressed as a fraction (between 0 and 1) of the x length "
                   "(after periodic repeats). The notation is: "
                   "start,end,increment or displacement_value.")
@click.option('-dy', '--displacement_y', default='0',
              help="Displacement to be made in the y (a2) direction. "
                   "See displacement_x for syntax")
@click.option('-spo', '--special_pointsonly', is_flag=True,
              help="Create only undistored and minimum energy displacements. "
                   "Overides the displacement_{x,y} commands. "
                   "May not be defined for all lattice/surface combinations. "
                   "Useful for benchmarking.")
@click.option('-prm', '--primitive', is_flag=True,
              help="Create primitive cell version (i.e non-orthogonal). ")
@click.option('-se', '--solute_elements', required=False,
              help="Solute element, created at unique distances from the stacking fault."
              " Will force the construction of a stable stacking fault."
              " Creates solutes at every unique layer away from the SF."
              " Can pass a list of elements using comma seperation"
              " E.g. Mg,Si,Cu. Can specify the creation of a vacancy using 'Vac'"
              " NOTE: will not generate symmetrically equivalent structures."
              " E.g. if Mg-Si solute solutes have been generated the script will skip Si-Mg")
@click.option('-msl', '--maxsolute_layer', default=None,
              help="Maximum layer to place solutes away from the SF")
@click.option('-tsl', '--testsolute_layer', is_flag=True,
              help="Place one solute at the midpoint (test) of the SF")
@click.option('-rsl', '--refsolute', is_flag=True,
              help="Place one solute at the origin of an undistorted slab with size of a SF")
@click.option('-sg', '--structure_group_name', required=True,
              help="Output AiiDA group to store created structures")
@click.option('-sgd', '--structure_group_description', default="",
              help="Description for output AiiDA group")
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints structures and extras but does not store anything")
def launch(lattice_size, matrix_element, lattice_and_surface,
           custom_structure,
           periodic_xrepeats, periodic_yrepeats, periodic_zrepeats,
           displacement_x, displacement_y, special_pointsonly,
           primitive, solute_elements, maxsolute_layer, testsolute_layer,
           refsolute, structure_group_name, structure_group_description,
           dryrun):
    """
    Script for creating stacking fault structures for a given size and matrix element. Generates
    a set of distorted structures using the 'tilted cell method', i.e. by adding fractional
    increments of the 'x' and 'y' cell vectors to the the 'z', vector.
    """
    STABLE_STACKING_NAME = 'stable_stacking'
    if not dryrun:
        structure_group = Group.get_or_create(
                             name=structure_group_name, description=structure_group_description)[0]
    else:
        structure_group = None



    extras = {}
    if lattice_and_surface:
        if lattice_size is None:
           raise Exception("Must specifiy a lattice_size if using lattice_and_surface")
        if matrix_element is None:
           raise Exception("Must specifiy a matrix_element if using lattice_and_surface")
        lattice_size = float(lattice_size)
        lattice_type, surface_plane = lattice_and_surface.split('_')
        surface_plane = "{"+str(surface_plane)+"}"
        orthogonal = not primitive
        extras = {
            'lattice_size':lattice_size,
            'lattice_type':lattice_type,
            'surface_plane':surface_plane,
            'matrix_element':matrix_element,
                      }
    special_points = {'undistorted':[0,0]}
    if lattice_and_surface == "FCC_111":
       xrepeats = periodic_xrepeats
       zrepeats = 3*periodic_zrepeats
       extras['z_direction']  = '<111>'
       if orthogonal:
           yrepeats = 2*periodic_yrepeats
           extras['x_direction']  = '<112>'
           extras['y_direction']  = '<110>'
           extras['orthogonal'] = 'True'
           special_points[STABLE_STACKING_NAME] = [0, 2./3.]
       else:
           yrepeats = periodic_yrepeats
           extras['x_direction']  = '<110>'
           extras['y_direction']  = '<110>'
           extras['orthogonal'] = 'False'
           special_points[STABLE_STACKING_NAME] = [1./3., 1./3.]
       undistorted_structure = ase.build.fcc111(
                                          matrix_element,
                                          [xrepeats,yrepeats,zrepeats],
                                          orthogonal=orthogonal,
                                          a=lattice_size)
    elif custom_structure:
        custom_structure = load_node(custom_structure)
        undistorted_structure = custom_structure.get_ase()
        extras = custom_structure.get_extras()
        if '_aiida_hash' in extras:
            del extras['_aiida_hash']
        #Check that x_direciton &  y_direction are already specified
        if 'label' not in extras:
            raise Exception("label not found in {} extras".format(custom_structure))
        if 'x_direction' not in extras:
            raise Exception("x_direction not found in {} extras".format(custom_structure))
        if 'y_direction' not in extras:
            raise Exception("y_direction not found in {} extras".format(custom_structure))
        if 'surface_plane' not in extras:
            raise Exception("surface_plane not found in {} extras".format(custom_structure))
        undistorted_structure =  undistorted_structure.repeat(
                                   [periodic_xrepeats, periodic_yrepeats, periodic_zrepeats]
                                                             )
    else:
       raise Exception("Could not process lattice_and_surface: {}".format(lattice_and_surface))


    #Add extras common to specified & custom lattice
    extras['periodic_xrepeats'] = periodic_xrepeats
    extras['periodic_yrepeats'] = periodic_yrepeats
    extras['periodic_zrepeats'] = periodic_zrepeats


    undistorted_structure.pbc = [True, True, True] # DFT structures always periodic
    a1 = undistorted_structure.get_cell()[0]/float(periodic_xrepeats)
    a2 = undistorted_structure.get_cell()[1]/float(periodic_yrepeats)
    a3 = undistorted_structure.get_cell()[2]/float(periodic_zrepeats)

    dispx_array = get_displacements_array(displacement_x)
    dispy_array = get_displacements_array(displacement_y)
    displacements = [[d_x, d_y] for d_x in dispx_array for d_y in dispy_array]
    special_pointnames = []

    if special_pointsonly:
        displacements = [] # overide any user displacements
        for sp_name in special_points:
            d_x, d_y = special_points[sp_name]
            displacements.append([d_x, d_y])
            special_pointnames.append(sp_name)

    if solute_elements:
        displacements = [] # overide any user displacements
        if STABLE_STACKING_NAME not in special_points:
            raise Exception("{} has no stable_stacking structure defined "
                            "".format(lattice_and_surface))
        if refsolute:
            d_x, d_y = special_points['undistorted']
            special_pointnames.append('undistorted')
        else:
            d_x, d_y = special_points[STABLE_STACKING_NAME]
            special_pointnames.append(STABLE_STACKING_NAME)
        displacements.append([d_x, d_y])

    for displacement in displacements:
        d_x, d_y = displacement
        extras['displacement_x'] = d_x
        extras['displacement_y'] = d_y
        if special_pointsonly:
            extras['special_point'] = special_pointnames.pop(0)
        distorted_structure = undistorted_structure.copy()
        distorted_structure.cell[2] += a1*d_x
        distorted_structure.cell[2] += a2*d_y
        store_asestructure(distorted_structure, extras, structure_group, dryrun)

    solute_elements = prep_elementlist(solute_elements)
    for solute_element in solute_elements:
        extras['sol1_element'] = solute_element
        layer_frame = get_layer_frame(distorted_structure, (0,0,1))
        layer_frame = layer_frame.drop_duplicates("layer_index").reset_index()
        solute_layers = range(int(len(layer_frame)/2))
        if refsolute:
            solute_layers = [0]
        if testsolute_layer:
            solute_layers = [int(len(layer_frame)/2)-1]
        for i in solute_layers:
            solute_structure = distorted_structure.copy()
            solute_index = int(layer_frame.loc[i]['structure_index'])
            solute_structure[solute_index].symbol = solute_element
            extras['sol1_index'] = solute_index
            extras['sol1sf_distance'] = layer_frame.loc[i]['layer_distance']
            extras['sol1layer_index'] = int(layer_frame.loc[i]['layer_distance'])
            store_asestructure(solute_structure, extras, structure_group, dryrun)
            if maxsolute_layer and i >= int(maxsolute_layer):
                break


if __name__ == "__main__":
   launch()
