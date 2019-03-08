#!/usr/bin/env python
import aiida
aiida.try_load_dbenv()
from aiida.orm.group import Group
from aiida_create_solutesupercell_structures import *
import ase
import ase.build
import click
import numpy as np

def get_displacements_array(displacement):
    if len(displacement.split(',')) == 3:
       d_min,d_max,d_inc = displacement.split(',')
       if dmax > 1:
           print "WARNING: max displacement {} is larger than 1".format(displacement)
       displacements = np.arange(d_min,d_max,d_inc)
    else:
       displacements = np.array([float(displacement)])
    return displacements


@click.command()
@click.option('-a', '--lattice_size', required=True,
              help="lattice length (in Ang) to use")
@click.option('-me', '--matrix_element', required=True,
              help="element to be used as the matrix")
@click.option('-l_surf', '--lattice_and_surface',
              type=click.Choice(["FCC_111"]),
              help="lattice and surface to be used. "
              "FCC_111: <112>(x) <110>(y) <111>(z)", required=True)
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
@click.option('-sg', '--structure_group_name', required=True,
              help="Output AiiDA group to store created structures")
@click.option('-sgd', '--structure_group_description', default="",
              help="Description for output AiiDA group")
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints structures and extras but does not store anything")
def launch(lattice_size, matrix_element, lattice_and_surface,
           periodic_xrepeats, periodic_yrepeats, periodic_zrepeats,
           displacement_x, displacement_y, special_pointsonly,
           structure_group_name, structure_group_description,
           dryrun):
    """
    Script for creating stacking fault structures for a given size and matrix element. Generates
    a set of distorted structures using the 'tilted cell method', i.e. by adding fractional
    increments of the 'x' and 'y' cell vectors to the the 'z', vector.
    """
    if not dryrun:
        structure_group = Group.get_or_create(
                             name=structure_group_name, description=structure_group_description)[0]
    else:
        structure_group = None

    lattice_size = float(lattice_size)
    lattice_type, surface_plane = lattice_and_surface.split('_')
    surface_plane = "{"+str(surface_plane)+"}"

    extras = {
        'lattice_size':lattice_size,
        'lattice_type':lattice_type,
        'surface_plane':surface_plane,
        'matrix_element':matrix_element,
                  }

    special_points = {'undistorted':[0,0]}
    if lattice_and_surface == "FCC_111":
       undistorted_structure = ase.build.fcc111(
                                            matrix_element,
                                           [periodic_xrepeats,
                                          2*periodic_yrepeats,
                                          3*periodic_zrepeats],
                                          orthogonal=True,
                                          a=lattice_size)
       special_points['stable_stacking'] = [2./3., 0]
       extras['x_direction']  = '<112>'
       extras['y_direction']  = '<110>'
       extras['z_direction']  = '<111>'
    else:
       raise Exception("Could not process lattice_and_surface: {}".format(lattice_and_surface))

    undistorted_structure.pbc = [True, True, True] # DFT structures always periodic
    a1 = undistorted_structure.get_cell()[0]/float(periodic_xrepeats)
    a2 = undistorted_structure.get_cell()[1]/float(periodic_yrepeats)
    a3 = undistorted_structure.get_cell()[2]/float(periodic_zrepeats)

    dispx_array = get_displacements_array(displacement_x)
    dispy_array = get_displacements_array(displacement_y)
    displacements = [[d_x, d_y] for d_x in dispx_array for d_y in dispy_array]

    if special_pointsonly:
        displacements = [] # overide any user displacements
        special_pointnames = []
        for sp_name in special_points:
            d_x, d_y = special_points[sp_name]
            displacements.append([d_x, d_y])
            special_pointnames.append(sp_name)

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

if __name__ == "__main__":
   launch()
