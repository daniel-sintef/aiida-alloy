#!/usr/bin/env python
from __future__ import print_function
import aiida
aiida.load_dbenv()
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm import Group, WorkCalculation
from aiida.orm.data.structure import StructureData
import click
from ase import units
import aiida_utils
import sys

# show default values in click
orig_init = click.core.Option.__init__


def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True


click.core.Option.__init__ = new_init

# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-wg', '--work_group', default=None,
              type=str, help="verdi work group to be converted to runner format")
@click.option('-sg', '--structure_group', default=None,
              type=str, help="verdi structure group to be converted to runner format")
@click.option('-f', '--filename', required=False, default=False,
         type=str, help="filename for outputfile, default = work_group+\".input.data\"")
def createjob(work_group,structure_group,filename):
    ''' e.g.
    ./aiida_export_group_to_runner.py -wg kmc_1000K_4
    ./aiida_export_group_to_runner.py -wg Al6xxxDB_passingsubset
    ./aiida_export_group_to_runner.py -sg Al6xxxDB_structuregroup
    '''
    if work_group:
        input_group = work_group
    elif structure_group:
        input_group = structure_group
    else:
        raise Exception("Must select at least one group")

    qb = QueryBuilder()
    qb.append(Group, filters={'name': input_group}, tag='g')
    if work_group:
        print('work_group:', work_group)
        print('filename  :',filename)
        qb.append(WorkCalculation, tag='job', member_of='g')
    if structure_group:
        print('structure_group:', structure_group)
        qb.append(StructureData, tag='structure', member_of='g')
    all_nodes = [x[0] for x in qb.all()]

    def get_workcalc_runnerdata(worknode):
        #TODO enable multi-structure support
        try:
            ase_structure = worknode.out.output_structure.get_ase()
        except Exception:
            ase_structure = worknode.inp.structure.get_ase()

        try:
            structure_path = worknode.inp.structure.get_extras()['structure_path']
        except Exception:
            structure_path = ''

        energy = worknode.out.output_parameters.get_attrs()['energy']  # units?

        #TODO: this section splits for SCF and relax, should fix & merge
        print("worknode: ", worknode)
        try:
            # SCF
            forces = worknode.out.output_array.get_array('forces')  # units?
        except Exception:
            # Relax (probably)
            forces = worknode.out.CALL.out.CALL.out.output_trajectory.get_array('forces')
        forces = forces[-1]

        #TODO: this section splits for SCF and relax, should fix & merge
        try:
            # SCF
            path = worknode.out.CALL.out.retrieved.get_abs_path()
        except Exception:
            # Relax (probably)
            path = "path support only for SCF calcs"
        return ase_structure, energy, forces, worknode.uuid, path, structure_path

    def get_structure_runnerdata(structurenode):
        ase_structure = structurenode.get_ase()

        try:
            structure_path = structurenode.get_extras()['structure_path']
        except Exception:
            structure_path = ''

        return ase_structure, structurenode.uuid, structure_path

    # This is from units tool, but ase uses other conversion factors
    #angstrom_to_bohrradius = 1.8897261
    #eV_to_Hartree = 0.036749325
    #eV_per_angstrom_to_hartree_per_bohrradius = 0.019446905

    fileOut = open(input_group+".input.data", "w")

    angstrom_to_bohrradius = 1./units.Bohr
    eV_to_Hartree = 1/units.Hartree
    eV_per_angstrom_to_hartree_per_bohrradius = units.Bohr/units.Hartree

    aiida_utils.create_READMEtxt()
    if filename == False:
        fileOut = open("aiida_exported_group_"+input_group+".input.data", "w")
    else:
        fileOut = open(filename, "w")

    for idx, node in enumerate(all_nodes):
        try:
            if work_group:
               ase_structure, energy, forces, uuid, path, structure_path = get_workcalc_runnerdata(node)
               print(idx, "ene (eV)", energy, uuid, path)
            if structure_group:
               ase_structure, uuid, structure_path = get_structure_runnerdata(node)
        except AttributeError:
            print('this worknode has errors')
            continue

        fileOut.write("begin\ncomment uuid: {}\n".format(uuid))
        fileOut.write("comment structure_path: {}\n".format(structure_path))

        # write the cell
        cell = ase_structure.cell*angstrom_to_bohrradius
        for idx_cell, i in enumerate(cell):
            fileOut.write("lattice %.10f %.10f %.10f\n" %
                          (cell[idx_cell][0], cell[idx_cell][1], cell[idx_cell][2]))

        #  write the positions
        nr_of_atoms = ase_structure.positions.shape[0]
        for idx_pos in range(nr_of_atoms):
            atCor = ase_structure.positions[idx_pos]*angstrom_to_bohrradius
            element = ase_structure.get_chemical_symbols()[idx_pos]
            if work_group:
                atFor = forces[idx_pos]*eV_per_angstrom_to_hartree_per_bohrradius
                fileOut.write("atom   %.6f    %.6f   %.6f %s  0.0   0.0  %.10f  %.10f  %.10f\n" %
                              (atCor[0], atCor[1], atCor[2],
                               element,
                               atFor[0], atFor[1], atFor[2]))
            elif structure_group:
                fileOut.write("atom   %.6f    %.6f   %.6f %s \n" %
                              (atCor[0], atCor[1], atCor[2], element))
            else:
                raise Exception("Must have either work or structure group")

        if work_group:
            fileOut.write("energy %.15f\n" % (energy*eV_to_Hartree))
        fileOut.write("charge 0\nend\n")

    fileOut.close()
    return


if __name__ == "__main__":
    createjob()
