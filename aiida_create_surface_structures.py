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
       displacements = np.arange(
                           float(d_min),float(d_max),float(d_inc)
                           )
    else:
       displacements = np.array([float(displacement)])
    return displacements


@click.command()
@click.option('-a', '--lattice_size', required=True,
              help="lattice length (in Ang) to use")
@click.option('-me', '--matrix_element', required=True,
              help="element to be used as the matrix")
@click.option('-l_surf', '--lattice_and_surface',
              type=click.Choice(["FCC_111", "FCC_110", "FCC_100"]),
              help="lattice and surface to be used. "
              "FCC_111: <112>(x) <110>(y) <111>(z) "
              "FCC_110: <100>(x) <1-10>(y) <110>(z) "
              "FCC_100: <100>(x) <010>(y) <001>(z) "
              , required=True)
@click.option('-pxr', '--periodic_xrepeats', default=1, type=int,
              help="periodic repeats in the x (a1) direction")
@click.option('-pyr', '--periodic_yrepeats', default=1, type=int,
              help="periodic repeats in the y (a2) direction")
@click.option('-pzr', '--periodic_zrepeats', default=1, type=int,
              help="periodic repeats in the z (a3) direction")
@click.option('-vact', '--vacuum_thickness', default='0',
              help="The thickness of the vacuum in Ang"
                   "The notation is: "
                   "start,end,increment or displacement_value.")
@click.option('-sg', '--structure_group_name', required=True,
              help="Output AiiDA group to store created structures")
@click.option('-sgd', '--structure_group_description', default="",
              help="Description for output AiiDA group")
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints structures and extras but does not store anything")
def launch(lattice_size, matrix_element, lattice_and_surface,
           periodic_xrepeats, periodic_yrepeats, periodic_zrepeats, vacuum_thickness,
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

    lattice_size = float(lattice_size)
    lattice_type, surface_plane = lattice_and_surface.split('_')
    surface_plane = "{"+str(surface_plane)+"}"

    extras = {
        'lattice_size':lattice_size,
        'lattice_type':lattice_type,
        'surface_plane':surface_plane,
        'matrix_element':matrix_element,
        'periodic_xrepeats': periodic_xrepeats,
        'periodic_yrepeats': periodic_yrepeats,
        'periodic_zrepeats': periodic_zrepeats
                  }

    if lattice_and_surface == "FCC_111":
       undistorted_structure = ase.build.fcc111(
                                            matrix_element,
                                           [periodic_xrepeats,
                                          2*periodic_yrepeats,
                                          3*periodic_zrepeats],
                                          orthogonal=True,
                                          a=lattice_size)
       extras['x_direction']  = '<112>'
       extras['y_direction']  = '<110>'
       extras['z_direction']  = '<111>'
    elif lattice_and_surface == "FCC_110":
       undistorted_structure = ase.build.fcc110(
                                            matrix_element,
                                           [periodic_xrepeats,
                                            periodic_yrepeats,
                                          2*periodic_zrepeats],
                                          a=lattice_size)
       extras['x_direction']  = '<100>'
       extras['y_direction']  = '<1-10>'
       extras['z_direction']  = '<110>'
    elif lattice_and_surface == "FCC_100":
       undistorted_structure = ase.build.fcc100(
                                            matrix_element,
                                           [periodic_xrepeats,
                                            periodic_yrepeats,
                                          2*periodic_zrepeats],
                                          a=lattice_size)
       extras['x_direction']  = '<100>'
       extras['y_direction']  = '<010>'
       extras['z_direction']  = '<001>'
    else:
       raise Exception("Could not process lattice_and_surface: {}".format(lattice_and_surface))

    undistorted_structure.pbc = [True, True, True] # DFT structures always periodic
    original_z_length = undistorted_structure.cell[2][2]

    vacuum_thicknesses = get_displacements_array(vacuum_thickness)

    for vct_i in vacuum_thicknesses:
        extras['vacuum_thickness'] = vct_i
        distorted_structure = undistorted_structure.copy()
        distorted_structure.cell[2][2] += vct_i
        store_asestructure(distorted_structure, extras, structure_group, dryrun)

if __name__ == "__main__":
   launch()
