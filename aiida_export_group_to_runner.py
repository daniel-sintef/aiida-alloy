#!/usr/bin/env python
from __future__ import print_function
import aiida
aiida.try_load_dbenv()
from aiida.orm import Node
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm import Calculation, Group, WorkCalculation
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.trajectory import TrajectoryData
from aiida.orm.utils import load_node, WorkflowFactory
import click
from ase import units
import aiida_utils
import sys
import numpy as np

#Define unit conversions
ANGSTROM_TO_BOHRRADIUS = 1./units.Bohr
EV_PER_ANGSTROM_TO_HARTREE_PER_BOHRRADIUS = units.Bohr/units.Hartree
EV_TO_HARTREE = 1/units.Hartree

def get_allnodes_fromgroup(group_name):
    qb = QueryBuilder()
    qb.append(Group, filters={'name': group_name}, tag='g')
    qb.append(Node, tag='job', member_of='g')
    all_nodes = [x[0] for x in qb.all()]
    return all_nodes

def write_runner_commentline(fileout, uuid, extra_comments={}):
    fileout.write("begin\ncomment ")
    fileout.write("uuid: {} ".format(uuid))
    for label in extra_comments:
        fileout.write("{}: {} ".format(label, extra_comments[label]))
    fileout.write("\n")
    return

def write_runner_cell(fileout, cell):
    cell = cell * ANGSTROM_TO_BOHRRADIUS
    for idx_cell, i in enumerate(cell):
        fileout.write("lattice %.10f %.10f %.10f\n" %
                      (cell[idx_cell][0], cell[idx_cell][1], cell[idx_cell][2]))
    return

def write_runner_atomlines(fileout, atomiccoord_array, elements, atomicforce_array=None):
    """
    Assumes input units of eV and angstrom
    """
    if atomicforce_array is None:
        atomicforce_array = np.zeros(atomiccoord_array.shape)

    atomiccoord_array = atomiccoord_array * ANGSTROM_TO_BOHRRADIUS
    atomicforce_array = atomicforce_array * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHRRADIUS

    for i in range(len(atomiccoord_array)):
        atomiccoord = atomiccoord_array[i]
        element = elements[i]
        atomicforce = atomicforce_array[i]

        fileout.write("atom   %.6f    %.6f   %.6f "
                      "%s  0.0   0.0  "
                      "%.10f  %.10f  %.10f\n" %
                      (atomiccoord[0], atomiccoord[1], atomiccoord[2],
                       element,
                       atomicforce[0], atomicforce[1], atomicforce[2]))
    return

def write_runner_finalline(fileout, energy=0, charge=0):
    fileout.write("energy %.15f\n" % (energy*EV_TO_HARTREE))
    fileout.write("charge %.15f\nend\n" % charge)
    return


def write_structure_torunner(fileout, structure_node, extra_comments={}):
    # get structure path, if applicable
    ase_structure = structure_node.get_ase()

    cell = ase_structure.get_cell()
    positions = ase_structure.get_positions()
    elements = ase_structure.get_chemical_symbols()

    write_runner_commentline(fileout, structure_node.uuid, extra_comments=extra_comments)
    write_runner_cell(fileout, cell)
    write_runner_atomlines(fileout, positions, elements)
    write_runner_finalline(fileout)
    return

def get_timesorted_basenodes(relaxworknode):
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": relaxworknode.uuid}, tag="relaxworknode")
    q.append(WorkCalculation, output_of="relaxworknode",
             project=["id", "ctime",  "*"],  tag="calc")
    q.order_by({"calc": "ctime"})
    timesorted_scf = [x[2] for x in q.all()]
    return timesorted_scf

def get_timesorted_scfs(worknode, relax_worknode=False):
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": worknode.uuid}, tag="worknode")
    output_tag = "worknode"
    if relax_worknode:
        output_tag = "worknode2"
        q.append(WorkCalculation, tag=output_tag, output_of="worknode")
    q.append(Calculation, output_of=output_tag, project=["id", "ctime",  "*"],  tag="calc")
    q.order_by({"calc": "ctime"})
    timesorted_scf = [x[2] for x in q.all()]
    return timesorted_scf

def write_pwbase_torunner(fileout, pwbasenode, extra_comments={}):
    scf_node = get_timesorted_scfs(pwbasenode)[-1]

    ase_structure = scf_node.inp.structure.get_ase()
    cell = ase_structure.get_cell()
    positions = ase_structure.get_positions()
    elements = ase_structure.get_chemical_symbols()

    atomicforce_array = scf_node.out.output_array.get_array('forces')[-1]
    energy = scf_node.out.output_parameters.get_attr('energy')

    write_runner_commentline(fileout, pwbasenode.uuid, extra_comments=extra_comments)
    write_runner_cell(fileout, cell)
    write_runner_atomlines(fileout, positions, elements, atomicforce_array=atomicforce_array)
    write_runner_finalline(fileout, energy=energy)
    return

def get_timesorted_trajectories(relaxworkcalc):
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": relaxworkcalc.uuid}, tag="relaxworkcalc")
    q.append(WorkCalculation, tag="baseworkcalc", output_of="relaxworkcalc")
    q.append(Calculation, output_of="baseworkcalc", tag="calc")
    q.append(TrajectoryData, output_of="calc", project=["id", "ctime",  "*"], tag="traj")
    q.order_by({"traj": "ctime"})
    timesorted_trajectories = [x[2] for x in q.all()]
    return timesorted_trajectories

def get_arraysbyname_fromtrajectories(timesorted_trajectories, arrayname):
    timesorted_arrays = [x.get_array(arrayname) for x in timesorted_trajectories]
    return np.concatenate(timesorted_arrays)

def write_pwrelax_torunner(fileout, relax_node, extra_comments={}):
    trajectories = get_timesorted_trajectories(relax_node)

    timesorted_cells = get_arraysbyname_fromtrajectories(trajectories, 'cells')
    timesorted_positions = get_arraysbyname_fromtrajectories(trajectories, 'positions')
    timesorted_forces = get_arraysbyname_fromtrajectories(trajectories, 'forces')
    if len(timesorted_forces) == 1:
        energy = relax_node.out.CALL.out.CALL.out.output_parameters.get_attr('energy')
        timesorted_energy = [energy]
    else:
        timesorted_energy = get_arraysbyname_fromtrajectories(trajectories, 'energy')
    elements = trajectories[0].get_array('symbols') # assume unchangin

    extra_comments={"trajectory_step":None}
    for i in range(len(timesorted_cells)):
        extra_comments["trajectory_step"] = i
        write_runner_commentline(fileout, relax_node.uuid, extra_comments=extra_comments)
        write_runner_cell(fileout, timesorted_cells[i])
        write_runner_atomlines(fileout,
           timesorted_positions[i], elements, atomicforce_array=timesorted_forces[i])
        write_runner_finalline(fileout, energy=timesorted_energy[i])

    if bool(relax_node.inp.final_scf):
        final_basenode = get_timesorted_basenodes(relax_node)[-1]
        extra_comments["trajectory_step"] = "final_scf"
        extra_comments["parent_uuid"] = relax_node.uuid
        write_pwbase_torunner(fileout, final_basenode, extra_comments=extra_comments)

    return


# show default values in click
orig_init = click.core.Option.__init__
def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True
click.core.Option.__init__ = new_init
# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-gn', '--group_name', default=None,
              type=str, help="Group to export identified by name")
@click.option('-f', '--filename', required=False, default=False,
         type=str, help="filename for outputfile, default = work_group+\".input.data\"")
@click.option('-sreadme', '--supress_readme', is_flag=True,
         help="supresses the generation of a readme file")
@click.option('-v', '--verbose', is_flag=True,
         type=str, help="Enables verbosity")
def createjob(group_name, filename, supress_readme, verbose):
    ''' e.g.
    ./aiida_export_group_to_runner.py -gn Al6xxxDB_structuregroup
    '''

    all_nodes = get_allnodes_fromgroup(group_name)
    if not supress_readme:
        aiida_utils.create_READMEtxt()

    if filename == False:
        fileout = open("aiida_exported_group_"+group_name+".input.data", "w")
    else:
        fileout = open(filename, "w")

    if verbose:
        print('filename  :', filename)
        print('structure_group:', group_name)

    for node in all_nodes:
        if verbose:
            print("Writing node: {}".format(node.uuid))
        if isinstance(node, StructureData):
            write_structure_torunner(fileout, node)
        elif isinstance(node, WorkCalculation):
            process_label = node.get_attrs()['_process_label']
            if process_label == "PwBaseWorkChain":
                write_pwbase_torunner(fileout, node)
            elif process_label == "PwRelaxWorkChain":
                write_pwrelax_torunner(fileout, node)
            else:
                print("Could not identify node, skipping")
        else:
            print("Could not identify node, skipping")

    fileout.close()
    return


if __name__ == "__main__":
    createjob()
